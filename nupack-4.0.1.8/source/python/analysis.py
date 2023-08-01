import functools, warnings, itertools, collections, numpy, pandas, decimal, scipy.sparse

from .rebind import Tuple, List, Dict, Union, Callable, forward
from .core import SetSpec, Structure, Complex, Sparsity, PairMatrix, Strand, complexes_from_max_size, Tube, Local, SequenceList, standardize_alphabet
from .utility import match, str_repr, long_output, printable, NamedTuple, \
    BaseResult, Loadable, from_argument, check_instance, check_instances, create_name
from .model import Model, Ensemble
from .rotation import lowest_rotation
from .concentration import solve_complex_concentrations
from .thermo import StructureEnergy, Job, CostsJob, PFJob, PairsJob, SuboptJob, SampleJob, ComputeOptions, submit

################################################################################

class ComplexResult(NamedTuple):
    '''Output result for analysis of a single Complex, with a lot of optional fields'''
    model: Model
    # output of pfunc
    pfunc: decimal.Decimal = None
    free_energy: float = None
    # output of mfe
    mfe_stack: float = None
    # output of ensemble_size
    ensemble_size: int = None
    # output of pairs
    pairs: PairMatrix = None
    # output of costs
    costs: numpy.ndarray = None
    # output of subopt
    mfe: list = None
    subopt: list = None
    # output of sample
    sample: list = None

    def structure_mfe(self) -> float:
        return self.mfe[0].energy if self.mfe else float('inf')

    def __repr__(self):
        return str_repr(self)

    def __str__(self):
        return  '{%s}' % ', '.join('{}: {}'.format(k, v) for k, v
            in zip(self._fields, self) if v is not None)

################################################################################

class Options(NamedTuple):
    '''Analysis options which apply to every complex'''

    num_sample: int = 0
    energy_gap: float = 0
    sparsity_fraction: float = 1.0
    sparsity_threshold: float = 0.0
    cache_bytes: float = -1.0

    from_argument = classmethod(from_argument)

################################################################################

def complex_analysis(complexes, model, *, compute, options=None):
    '''
    Analyze a set of complexes for their complex ensemble thermodynamic properties
    - complexes: a ComplexSet, Tube, or list of complexes
    - compute: a list of physical quantities to compute [pfunc, mfe, pairs, sample, subopt, ensemble_size]
    - options: a dict or instance of `Options`: the options to apply to each calculation
    - model: the free energy model to use
    '''
    complexes = check_instances(set(complexes), Complex)
    jobs, counts = [], []
    options = Options.from_argument(options)

    # If the model is rnadna, standardize complexes
    if model.material == 'DNA/RNA':
        complexes_s = standardize_alphabet(complexes, model.alphabet())
        complexes_s = [complexes_s[i] for i in range(0,len(complexes_s))]
    else:
        complexes_s = complexes

    for c in complexes_s:
        for task in compute:
            if task == 'pfunc':
                jobs.append(Job(c, PFJob()))
            elif task == 'mfe':
                jobs.append(Job(c, SuboptJob(gap=0)))
            elif task == 'costs':
                jobs.append(Job(c, CostsJob()))
            elif task == 'pairs':
                jobs.append(Job(c, PairsJob(Sparsity(
                    row_size=int(round(options.sparsity_fraction * c.nt())),
                    threshold=options.sparsity_threshold))))
            elif task == 'sample':
                jobs.append(Job(c, SampleJob(number=options.num_sample)))
            elif task == 'subopt':
                jobs.append(Job(c, SuboptJob(gap=options.energy_gap)))
            elif task == 'ensemble_size':
                counts.append(Job(c, PFJob()))
            else:
                raise KeyError('%r is not a valid analysis task' % task)

    compute_ops = ComputeOptions(max_bytes=options.cache_bytes)

    if 'ensemble_size' in compute:
        cmodel = model.copy()
        cmodel.beta = 0
    else:
        cmodel = None
    results = [j and submit(j, model=m, options=compute_ops) for j, m in zip([jobs, counts], [model, cmodel])]
    results, counts = [r and r.get() for r in results]

    output = {}

    if 'sample' in compute:
        samples = {}
        for k, v in results.items():
            nicks = k.nicks()
            if v.sample:
                samples[k] = [Structure(s, nicks) for s in v.sample.value().structures]

    has_non_count = bool(set(compute).difference(['ensemble_size']))

    for i in range(0, len(complexes)):
        c = complexes_s[i]
        c_orig = complexes[i]
        key = SequenceList(c)
        v = results[key] if has_non_count else None

        sym = c.symmetry()
        bonus = c.bonus
        # now have to fix up all this distinguishability ... stuff
        r = {}

        if 'mfe' in compute or 'subopt' in compute or 'costs' in compute:
            assert v.subopt, c
            nicks = c.nicks()
            strucs = v.subopt.take()
            strucs = [StructureEnergy(structure=Structure(s.structure, nicks),
                stack_energy=s.stack_energy+bonus, energy=s.energy + bonus) for s in strucs.structures]
            r['mfe_stack'] = v.mfe.take().energy + bonus

            mfe = min(s.energy for s in strucs) if strucs else None
            r['mfe'] = [s for s in strucs if s.energy <= mfe + 1e-4]

            if 'subopt' in compute:
                r['subopt'] = strucs
            
            if 'costs' in compute:
                r['costs'] = v.costs.take().matrix.cast(numpy.ndarray)

        if 'ensemble_size' in compute:
            pf = counts[key].pfunc
            assert pf, c
            r['ensemble_size'] = int(round(decimal.Decimal(pf.take().logq).exp())) // sym

        if 'pfunc' in compute or 'sample' in compute or 'pairs' in compute:
            assert v.pfunc, c
            pf = v.pfunc.take()
            factor = numpy.exp(-model.beta * bonus)
            r['pfunc'] = decimal.Decimal(pf.logq).exp() * decimal.Decimal(factor / sym)
            r['free_energy'] = (numpy.log(float(sym)) - pf.logq) / model.beta

        if 'sample' in compute:
            x = samples[key]
            r['sample'] = x[:options.num_sample]
            samples[key] = x[options.num_sample:]

        if 'pairs' in compute:
            assert v.pairs, c
            r['pairs'] = v.pairs.take().matrix

        output[c_orig] = ComplexResult(model=model, **r)

    return Result(complexes=output)

################################################################################

def into_blocks(matrix, lengths):
    '''split a matrix into blocks with respect to the given lengths'''
    prefix = numpy.cumsum(lengths)
    return [[matrix[p1-l1:p1, p2-l2:p2] for p2, l2 in zip(prefix, lengths)]
                                        for p1, l1 in zip(prefix, lengths)]

################################################################################

def ensemble_pair_fractions(strands, complex_data):
    '''
    Calculate ensemble pair fractions given the pair probability matrices for a set of complexes
    - strands: list of Strand
    - complex_data: map from Complex to pair probability matrix (ndarray or sparse)
    Should work with either sparse or non-sparse matrices
    '''
    strands = tuple(strands)
    assert len(set(strands)) == len(strands), 'Strands should be unique'
    P = [[0 for _ in strands] for _ in strands]
    for k, (c, x) in complex_data.items():
        blocks = into_blocks(x, tuple(map(len, k)))
        for i, row in zip(k, blocks):
            for j, p in zip(k, row):
                P[strands.index(i)][strands.index(j)] += c * p
    P = scipy.sparse.vstack(list(map(scipy.sparse.hstack, P)), format='csc')
    P.data /= numpy.asarray(P.sum(1)).reshape(-1)[P.indices]
    return P

################################################################################

def anonymous_complex(strands, alphabet=None):
    if alphabet is not None:
        strands = [Strand(s, name=str(i), alphabet=alphabet) for i, s in enumerate(SequenceList(strands, alphabet=alphabet))]
    else:
        strands = [Strand(s, name=str(i)) for i, s in enumerate(SequenceList(strands))]
    return Complex(strands)

################################################################################

def pfunc(strands, model) -> Tuple[float, float]:
    '''
    Calculate partition function of a single complex
    - model(Model): free energy model to use
    - Returns a tuple of (pfunc, free_energy)
    '''
    strands = anonymous_complex(strands, alphabet = model.alphabet())
    res = complex_analysis([strands], model, compute=['pfunc'])[strands]
    return (res.pfunc, res.free_energy)

################################################################################

def pairs(strands, model, *, sparsity_fraction=1, sparsity_threshold=0) -> PairMatrix:
    '''
    Calculate equilibrium pair probabilities of a single complex
    - model(Model): free energy model to use
    - sparsity_fraction: maximum fraction of each row to keep
    - sparsity_threshold: minimum pair probability to preserve
    - Returns a PairMatrix object
    '''
    strands = anonymous_complex(strands, alphabet = model.alphabet())
    res = complex_analysis([strands], model, compute=['pairs'],
        options=dict(sparsity_fraction=sparsity_fraction, sparsity_threshold=sparsity_threshold))
    return res[strands].pairs

################################################################################

def fraction_bases_unpaired(strands, model, max_size) -> float:
    '''
    Compute equilibrium fraction of unpaired bases
    - strands(dict): dict from strand (e.g. str) to its concentration
    - model(Model): free energy model to use
    '''
    tube = Tube({Strand(k, name=str(i)): v for i, (k, v) in enumerate(strands.items())}, 
        complexes=SetSpec(max_size), name='t')
    res = tube_analysis([tube], model=model, compute=['pairs'], options={'sparsity_fraction': 0})
    return res[tube].fraction_bases_unpaired

################################################################################

def mfe(strands, model) -> List[StructureEnergy]:
    '''
    Calculate MFE structures and energies for a single complex
    - model(Model): free energy model to use
    - Returns a list of StructureEnergy objects
    '''
    strands = anonymous_complex(strands, alphabet = model.alphabet())
    return complex_analysis([strands], model, compute=['mfe'])[strands].mfe

################################################################################

def energy(strands, structure, model) -> float:
    '''
    Calculate energy of a single structure of a complex
    - model(Model): free energy model to use
    '''
    bonus = getattr(strands, 'energy', 0.0)
    strands = anonymous_complex(strands, alphabet = model.alphabet())
    return model.structure_energy(strands, structure) + bonus

################################################################################

def structure_probability(strands, structure, model) -> float:
    '''
    Calculate equilibrium probability of a structure for a given complex
    - model(Model): free energy model to use
    '''
    strands = anonymous_complex(strands, alphabet = model.alphabet())
    res = complex_analysis([strands], model, compute=['pfunc'])[strands]
    e = energy(strands, structure, model)
    return numpy.exp(-model.beta * (e - res.free_energy))

prob = structure_probability

################################################################################

def ensemble_size(strands, model):
    '''
    Calculate number of secondary structures in a given complex ensemble
    - model(Model): free energy model to use
    '''
    strands = anonymous_complex(strands, alphabet = model.alphabet())
    res = complex_analysis([strands], model, compute=['ensemble_size'])
    return res[strands].ensemble_size

################################################################################

def subopt(strands, energy_gap, model):
    '''
    Calculate suboptimal structures falling below a specified energy gap of the MFE
    - model(Model): free energy model to use
    - energy_gap: suboptimal structure energy gap in kcal/mol
    '''
    strands = anonymous_complex(strands, alphabet = model.alphabet())
    res = complex_analysis([strands], model, compute=['subopt'], options={'energy_gap': energy_gap})
    return res[strands].subopt

################################################################################

def sample(strands, num_sample, model):
    '''
    Calculate a set of num_sample Boltzmann sampled structures
    - model(Model): free energy model to use
    - num_sample(int): number of samples
    '''
    strands = anonymous_complex(strands, alphabet = model.alphabet())
    res = complex_analysis([strands], model, compute=['sample'], options={'num_sample': num_sample})
    return res[strands].sample

################################################################################

class ConcentrationResult:
    '''
    Simple wrapper fo the equilibrium complex concentrations
    Use the `complex_concentrations` property to get a dict of these concentrations
    '''
    def __init__(self, strands, complexes, concentrations):
        self.strands = tuple(strands)
        self.complexes = tuple(complexes)
        self.concentrations = numpy.asarray(concentrations)
        assert len(self.complexes) == len(self.concentrations)

    @property
    def complex_concentrations(self):
        return dict(zip(self.complexes, self.concentrations))

################################################################################

class ConcentrationSolver:
    def __init__(self, strands, complex_results, *, distinguishable):
        '''
          Initialize the solver from
          - strands: an ordered list of Strand
        - - complex_results: a dict of Complex to ComplexResult
          Note that each ComplexResult must contain the partition function
          If distinguishable=True, then a Strand entered twice
          is equivalent to saying that the first entry is labeled differently than
          the second.
        '''
        self.strands = tuple(strands)
        self.temperature = tuple(complex_results.values())[0].model.temperature

        self.complexes = {k : float(v.pfunc.ln())
            for k, v in complex_results.items() if v.pfunc is not None}

        self.distinguishable = bool(distinguishable)

    def compute(self, concentrations, complexes=None, as_strands=True, **kws):
        '''
        concentrations: array
            list of molarity of each strand species
        complexes: None, int, or List[List[str]]
            None: calculate concentrations in an ensemble containing all prior calculated complexes
            int: calculate concentration in an ensemble containing complexes of only up to this size
            List[List[str]]: calculate concentrations in an ensemble of only these complexes
        as_strands: bool
            Whether initial concentrations are given for each strand or for each complex
        **kws:
            Custom options for the concentration solving algorithm (see solve_complex_concentrations)
        '''
        complexes = tuple(self.complexes.keys()) if complexes is None else complexes

        if as_strands:
            if len(concentrations) != len(self.strands):
                raise ValueError('strand number %d != concentration number %d' % (len(self.strands), len(concentrations)))
        else:
            if len(concentrations) != len(complexes):
                raise ValueError('complex number %d != concentration number %d' % (len(complexes), len(concentrations)))

        logq = [self.complexes[k] for k in complexes]
        indices = [[self.strands.index(s) for s in k] for k in complexes]
        x = solve_complex_concentrations(indices, logq, concentrations,
            kelvin=self.temperature, as_strands=as_strands,
            rotational_correction=self.distinguishable, **kws)
        return ConcentrationResult(self.strands, complexes, x)

################################################################################

@printable
class Result(BaseResult, Loadable):
    '''Class holding an analysis result'''

    names = ('tubes', 'complexes')

    def __init__(self, *, tubes=None, complexes=None):
        self.tubes = {} if tubes is None else dict(tubes)
        self.complexes = {} if complexes is None else dict(complexes)

    @property
    def fields(self):
        return (self.tubes, self.complexes)

    _unicode_cols = ['Complex', 'Pfunc', '\u0394G (kcal/mol)', 'MFE (kcal/mol)', 'Ensemble size']
    _pretty_cols = ['Complex', 'Pfunc', 'dG (kcal/mol)', 'MFE (kcal/mol)', 'Ensemble size']
    _raw_cols = ['complex', 'pfunc', 'free_energy', 'mfe', 'ensemble_size']

    _formats = [
        lambda c: c.name, # complex
        '{:.4e}'.format, # pfunc
        '{:.3f}'.format, # free energy
        '{:.3f}'.format, # mfe
        lambda i: str(i) if i < 1e10 else '{:.3e}'.format(i), # count
    ]

    _tube_fmts = [lambda x: '' if numpy.isnan(x) else '{:.3e}'.format(x)]

    def complex_frame(self, *, include_null=True, pretty: bool, unicode=True):
        '''pandas DataFrame with complexes'''
        data = [[k, v.pfunc, v.free_energy, v.mfe[0].energy if v.mfe else None, v.ensemble_size]
            for k, v in self.complexes.items()]
        if pretty:
            cols = self._unicode_cols if unicode else self._pretty_cols
        else:
            cols = self._raw_cols
        df = pandas.DataFrame(data, columns=cols)

        if include_null:
            return df

        df = df[[c for c in df.columns if not df[c].isnull().all()]]
        # sort by #strands, then alphabetically but move ( to the back
        keys = [(len(c), c.name.replace('(', '~')) for c in df[df.columns[0]]]
        return df.reindex(sorted(range(len(keys)), key=keys.__getitem__)).reset_index(drop=True)

    def tube_frame(self, *, pretty: bool):
        '''
        Return a pandas DataFrame with tubes
        - pretty: use spaces and capitalization
        '''
        if len(set(t.name for t in self.tubes)) != len(self.tubes):
            raise ValueError('Tube names are not unique')

        if self.complexes:
            complexes = self.complexes
        else:
            complexes = set(x for t in self.tubes for x in t.complexes)

        dfs = []
        for t, v in self.tubes.items():
            if pretty:
                df = pandas.DataFrame([[*x, numpy.nan] for x in v.complex_concentrations.items()],
                    columns=['Complex', t.name + ' (M)', ' '] )
            else:
                df = pandas.DataFrame(list(v.complex_concentrations.items()),
                    columns=['complex', t.name])
            # Sort by concentrations
            dfs.append(df.sort_values(by=df.columns[1], ascending=False).reset_index(drop=True))
        return pandas.concat(dfs, axis=1)

    def __str__(self):
        items = []
        with long_output():
            if self.complexes:
                df = self.complex_frame(include_null=False, pretty=True, unicode=False)
                fmts = dict(zip(df.columns, self._formats))
                string = df.to_string(formatters=fmts, na_rep='')
                items.append('Complex results:\n' + string)
            if any(len(t.complexes) > 1 for t in self.tubes):
                df = self.tube_frame(pretty=True)
                f = {k: self._tube_fmts[0] for k in df.columns}
                f[df.columns[0]] = lambda c: c.name
                items.append('Concentration results:\n' + df.to_string(formatters=f, na_rep='', index=False))
        return '\n'.join(items)

    def _repr_html_(self):
        '''Return html string for output in notebooks'''
        items = []
        with long_output():
            if self.complexes:
                df = self.complex_frame(include_null=False, pretty=True, unicode=True)
                f = dict(zip(df.columns, self._formats))
                string = df.to_html(formatters=f, na_rep='', index=False)
                items.append('<b>Complex results:</b> ' + string)
            if any(len(t.complexes) > 1 for t in self.tubes):
                df = self.tube_frame(pretty=True)
                f = {k: self._tube_fmts[0] for k in df.columns}
                f[df.columns[0]] = lambda c: c.name
                string = df.to_html(formatters=f, na_rep='', index=False)
                items.append('<b>Concentration results:</b> ' + string)
        return ' '.join(items)

################################################################################

class PairFractions:
    '''
    Simple wrapper class for ensemble pair fractions
    Use `to_array` or `to_sparse` to convert to numpy/scipy arrays.
    '''
    def __init__(self, sparse_matrix):
        self.matrix = sparse_matrix

    def to_array(self):
        '''Return numpy ndarray of the pair fractions'''
        return self.matrix.todense()

    def to_sparse(self):
        '''Return sparse array of the pair fractions'''
        return self.matrix

    def __str__(self):
        return str(self.to_array())

################################################################################

class TubeResult:
    '''
    Contains three fields:
    - complex_concentrations: a dict from complex to its equilibrium concentration
    - ensemble_pair_fractions: optionally, the equilibrium strand pairing matrix
    - fraction_bases_unpaired: optionally, the fraction of unpaired nucleotides in the test tube at equilibrium
    The latter two fields will only be populated if pair probability information has been computed.
    '''

    def __init__(self, complex_concentrations, *, tube, result=None):
        self.complex_concentrations = dict(complex_concentrations)

        if result is None or any(result[x].pairs is None for x in tube.complexes):
            self.ensemble_pair_fractions = None
            self.fraction_bases_unpaired = None
        else:
            self.fraction_bases_unpaired = sum(result[x].pairs.diagonal.sum() * complex_concentrations[x] for x in tube.complexes) \
                                   / sum(len(result[x].pairs.diagonal)  * complex_concentrations[x] for x in tube.complexes)
            raw = ensemble_pair_fractions(tube.strands,
                {x: (complex_concentrations[x], result[x].pairs.to_sparse()) for x in tube.complexes})
            self.ensemble_pair_fractions = PairFractions(raw)

    def __getitem__(self, complex):
        '''Get complex concentration by complex or complex name'''
        if not isinstance(complex, str):
            return self.complex_concentrations[complex]
        try:
            return next(v for k, v
                in self.complex_concentrations.items() if k == complex)
        except StopIteration:
            raise KeyError(complex) from None

    def __str__(self):
        return str(self.complex_concentrations)

################################################################################

def tube_analysis(tubes, model, *, compute=(), **kws):
    '''
    Analyze a set of tubes for their complex ensemble and concentration properties
    '''
    compute = ['pfunc'] + list(compute)
    result = complex_analysis(Tube.union(tubes), compute=compute, model=model, **kws)
    for t in tubes:
        solver = ConcentrationSolver(t.strands,
            {k: result[k] for k in t.complexes}, distinguishable=False)
        solved = solver.compute(t.concentrations)
        result.tubes[t] = TubeResult(solved.complex_concentrations, tube=t, result=result)
    return result

################################################################################

def complex_concentrations(tube, data, concentrations=None):
    '''Calculate the equilibrium concentrations for a given tube using already calculated complex results'''
    if concentrations is None:
        conc = tube.concentrations
    else:
        conc = [concentrations[strand] for strand in tube.strands]
    name = getattr(tube, 'name', 'tube')
    tube = Tube(zip(tube.strands, conc), complexes=tube.complexes, name=name)
    solver = ConcentrationSolver(tube.strands,
        {k: data[k] for k in tube.complexes}, distinguishable=False)

    res = solver.compute(tube.concentrations)
    return Result(tubes={tube: TubeResult(res.complex_concentrations, tube=tube)})

################################################################################
