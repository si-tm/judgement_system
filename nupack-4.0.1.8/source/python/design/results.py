from typing import Tuple, List, Dict, Union, Callable, NamedTuple
from scipy.sparse import csc_matrix
from decimal import Decimal
import pandas, numpy

from ..core import JSON, Pickleable, Sequence, Domain, Strand, Complex, SetSpec, Tube
from ..utility import printable, BaseResult, Loadable, long_output
from ..analysis import Result as AnalysisResult, ComplexResult as AnalysisComplexResult, TubeResult as AnalysisTubeResult
from .weights import style

################################################################################

class Partition:
    """Partition of complexes in design between active and passive sets"""
    mask: Tuple[bool, ...]
    deflate: float

################################################################################

class Stats:
    """Run statistics of a finished design"""
    num_leaf_evaluations: int
    num_reseeds: int
    num_redecompositions: Tuple[int, ...]
    offtargets_added_per_refocus: Tuple[int, ...]
    design_time: float
    analysis_time: float
    final_Psi: Partition
    seed: int

################################################################################

class ComplexResult:
    """Result of a complex calculation"""
    name: str
    # sequence: 'core.Complex'
    # structure: 'components.Structure'
    log_pfunc: float
    pair_probabilities: csc_matrix
    defect: float
    normalized_defect: float

################################################################################

class TubeComplex:
    """Result of a complex in a tube"""
    name: str
    concentration: float
    target_concentration: float
    nucleotide_defect: float
    structural_defect: float
    concentration_defect: float
    defect: float

################################################################################

class TubeResult:
    """Result of a tube calculation"""
    name: str
    nucleotide_concentration: float
    nucleotide_defect: float
    normalized_defect: float
    defect: float
    # complexes: 'std.Vector' # Tuple['TubeComplex', ...]

################################################################################

class Single:
    """One of the pareto optimal results at end of design"""
    domains: Dict[str, Sequence]
    strands: Dict[str, Sequence]
    defects: List[float]
    weighted_defects: List[float]

################################################################################

class RawResult(Pickleable):
    """Result of a design"""
    # model: 'components.ModelSettings'
    # parameters: 'components.Parameters'
    # stats: Stats'
    # complexes: 'std.Vector' # Tuple['Complex', ...]
    # tubes: 'std.Vector' # 'std.Vector' # Tuple['Tube', ...]
    success: bool

    def __init__(self, json=None, _fun_=None):
        """
        create a new Result, loading from JSON if specified

        Arguments:
            json: a valid nupack.JSON object specifying a Result object.
        """
        if json:
            self.move_from(self.new_from_json(json))
        else:
            _fun_(self)

################################################################################

@printable
class Mapping(BaseResult):
    '''Class for mapping from undesigned to designed entities'''

    names = ('tubes', 'complexes', 'strands', 'domains')

    def __init__(self):
        self.tubes = {}
        self.complexes = {}
        self.strands = {}
        self.domains = {}

    @property
    def fields(self):
        return (self.tubes, self.complexes, self.strands, self.domains)

    def domain_table(self, complements=True):
        return pandas.DataFrame([(k.name, str(v))
            for k, v in self.domains.items() if complements or not k.name.endswith('*')],
            columns=('Domain', 'Sequence'))

    def rules(self):
        return self.domains

    def __call__(self, key):
        for x in self.fields:
            y = x.get(key)
            if y is not None:
                return y
        raise KeyError(key)

    def get(self, key, default=None):
        '''Lookup a given domain'''
        return self.domains.get(key, default)

    def strand_table(self):
        return pandas.DataFrame([(k.name, str(v)) for k, v in self.strands.items()],
            columns=('Strand', 'Sequence'))

    def _fmt(self, html, x):
        if html:
            h = x.to_html(formatters={'Sequence': '1@@@@{}2@@@@'.format}, index=False)
            return h.replace('1@@@@', '<pre>').replace('2@@@@', '</pre>')
        else:
            return x.to_string(index=False)

    def __str__(self):
        with long_output():
            return 'Domain results:\n{}\n\nStrand results:\n{}'.format(
                self._fmt(False, self.domain_table()),
                self._fmt(False, self.strand_table()))

    def _repr_html_(self):
        with long_output():
            return '<b>Domain results</b>:\n{}<b>Strand results</b>:\n{}'.format(
                self._fmt(True, self.domain_table()),
                self._fmt(True, self.strand_table()))

################################################################################

@printable
class ComplexDefect(NamedTuple):
    '''Information on a complex's defect independent of its concentration'''
    defect: float
    normalized: float

################################################################################

@printable
class TubeDefect(NamedTuple):
    '''Information on a tube's defect'''
    defect: float
    normalized: float
    complexes: dict

################################################################################

@printable
class TargetDefect(NamedTuple):
    defect: float
    normalized: float
    concentration: float
    structural: float

################################################################################

@printable
class Defects:
    '''Complete breakdown of a design's resultant defect by tube and complex'''

    def __init__(self, defects, weighted, complexes, tubes, success, soft_constraints):
        objs = list(zip(['ensemble_defect'] + [c.soft_kind for c in soft_constraints],
            weighted, defects, ['Weighted ensemble defect'] + [c.soft_name for c in soft_constraints]))
        self.objectives = pandas.DataFrame(objs, columns=['kind', 'weighted', 'defect', 'name'])
        self.success = bool(success)

        self.complexes = pandas.DataFrame([(c.name, v.defect, v.normalized, c)
            for c, v in complexes.items() if len(c.structure)], columns=['complex_name', 'defect', 'normalized', 'complex'])

        self.tubes = pandas.DataFrame([(t.name, v.defect, v.normalized, t) for t, v in tubes.items()],
            columns=['tube_name', 'defect', 'normalized', 'tube'])

        self.tube_complexes = pandas.DataFrame([
            (t.name, c.name, r.structural, r.concentration, r.defect, t, c)
            for t, v in tubes.items() for c, r in v.complexes.items() if c in t.on_targets],
            columns=['tube_name', 'complex_name', 'structural', 'concentration', 'total', 'tube', 'complex'])

    @property
    def ensemble_defect(self):
        return self.objectives.defect[0]

    @property
    def weighted_ensemble_defect(self):
        return self.objectives.weighted[0]

    def _is_complex_design(self):
        return numpy.all(self.tube_complexes.concentration == 0)

    def _formatted(self, html, lns):
        g = '{:#.3g}'.format
        f = '{:#.3g}'.format

        # s =  [lns[0], style(pandas.DataFrame([[self.ensemble_defect, self.weighted_ensemble_defect]]),
        #    ['Ensemble defect', 'Weighted ensemble defect'], html)]
        # s = [lns[0], g(self.ensemble_defect)]

        df = self.objectives.groupby('name').sum(numeric_only=True).reset_index()
        df = df.sort_values('name', key=lambda x : x != 'Weighted ensemble defect')
        if len(df) > 1:
            df.loc[len(df)] = ['Total', df.weighted.sum(), 0]
        s = [lns[0], style(df, ['Objective type', 'Value'], html, weighted=f)]
        s += [lns[1], f(self.ensemble_defect)]
        # if len(self.objectives) > 1:
            # df = pandas.DataFrame([['Objective function + soft constraints', self.objectives.weighted.sum()]])
            # s += [lns[2], style(df, ['Augmented objective function', 'Weighted ensemble defect + weighted penalties'], html)]

        s += [lns[2], style(self.complexes, ['Complex', 'Complex defect (nt)', 'Normalized complex defect'],
            html, defect=g, normalized=lambda x: '{:#.3g}'.format(float(x)))]

        if not self._is_complex_design():
            s += [lns[3], style(self.tubes, ['Tube', 'Tube defect (M)', 'Normalized tube defect'], html, defect=g, normalized=f)]
            s += [lns[4], style(self.tube_complexes,
                ['Tube', 'On-target complex', 'Structural defect (M)', 'Concentration defect (M)', 'Total defect (M)'],
                html, concentration=g, total=g, structural=g)]

        return ''.join(s)

    def __str__(self):
        with long_output():
            return self._formatted(False, [
                'Objective function:\n', '\n\nEnsemble defect: ', '\n\n', '\n\nOn-target complex defects:\n',
                '\n\nTube defects:\n', '\nComplex contributions to tube defects:\n'
            ])

    def _repr_html_(self):
        with long_output():
            return self._formatted(True, [
                '<b>Objective function</b>:', '<b>Ensemble defect</b>: ', '<br><br><b>On-target complex defects:</b>',
                '<b>Tube defects:</b>', '<b>Complex contributions to tube defects:</b>'
            ])

################################################################################

class TargetConcentration(NamedTuple):
    actual: float
    target: float

################################################################################

class Concentrations:
    def __init__(self, tubes):
        self.table = pandas.DataFrame([(t.name, c.name, r.actual, r.target, c.nt(), t, c)
            for t, v in tubes.items() for c, r in v.items()],
            columns=['tube_name', 'complex_name', 'concentration', 'target_concentration', 'nucleotides', 'tube', 'complex'])

    def _split(self, html):
        on = self.table[self.table.target_concentration != 0]
        df = self.table.copy()
        # df['mass'] = df.target_concentration * df.nucleotides
        df['thr'] = 1e-2 * df.groupby('tube').transform('max').concentration
        off = df[(df.concentration >= df.thr) & (df.target_concentration == 0)].copy()
        empty = '\u2014' if html else '-'
        for t in set(df.tube).difference(off.tube):
            off.loc[len(off)] = [t.name, empty, float('nan'), 0, 0, t, None, 0]
        f = '{:#.3g}'.format
        return (style(on, ['Tube', 'Complex', 'Concentration (M)', 'Target concentration (M)'], html, concentration=f, target_concentration=f),
                style(off, ['Tube', 'Complex', 'Concentration (M)'], html, concentration=lambda x: empty if numpy.isnan(x) else f(x)))

    def __str__(self):
        with long_output():
            return 'On-target complex concentrations:\n%s\n\nSignificant off-target complex concentrations (>= 1%% max complex concentration in tube):\n%s' % self._split(False)

    def _repr_html_(self):
        with long_output():
            return '<b>On-target complex concentrations</b>:\n%s\n\n<b>Significant off-target complex concentrations (\u2265 1%% max complex concentration in tube)</b>:\n%s' % self._split(True)

################################################################################

@printable
class Result(NamedTuple):
    to_analysis: Mapping
    defects: Defects
    concentrations: Concentrations
    analysis_result: AnalysisResult
    stats: dict
    raw: 'SingleResult'
    design: 'Design'

    @property
    def ensemble_defect(self):
        return self.defects.ensemble_defect

    @property
    def domains(self):
        return self.to_analysis.domains

    def evaluate_with(self, domains=None, tubes=None, model=None,
        soft_constraints=None, defect_weights=None, objective_weight=None):
        '''
        Evaluate the same design with a differently specified set of tubes,
        weights, soft_constraints, defect_weights, or objective_weight
        '''
        # Apply non-domain customizations
        new = self.design.apply(domains=self.domains.values(), tubes=tubes, model=model,
            soft_constraints=soft_constraints, defect_weights=defect_weights,
            objective_weight=objective_weight)
        # Apply manually chosen domains
        if domains is not None:
            new = new.apply(domains)
        return new.evaluate()

    def _present(self):
        if self.defects._is_complex_design():
            return self.to_analysis, self.defects
        else:
            return self.to_analysis, self.defects, self.concentrations

    def _repr_html_(self):
        return ''.join(x._repr_html_() for x in self._present())

    def __repr__(self):
        return '<DesignResult: defect=%f>' % self.ensemble_defect

    def __str__(self):
        return '\n\n'.join(map(str, self._present()))

    @classmethod
    def build(cls, spec, raw_result):
        '''Build a Result object from a results.RawResult object'''
        r = raw_result
        lu = spec.lookup()
        D, S, X, T = [lambda k, d=d: d[getattr(k, 'name', k)] for d in lu] 

        # Fix reported ensemble defect to always be unweighted
        unweighted = [
            sum(t.normalized_defect for t in r.results[0].tubes) / len(r.results[0].tubes),
            *r.results[0].defects[1:]
        ]

        # Build Defects object from the raw outputs
        defects = Defects(
            defects=unweighted,
            weighted=r.results[0].weighted_defects,
            complexes={X(x): ComplexDefect(defect=x.defect, normalized=x.normalized_defect)
                for x in r.results[0].complexes if x.name in lu.complexes},
            tubes={T(t): TubeDefect(float(t.defect), t.normalized_defect,
                {X(x): TargetDefect(float(x.defect), float(x.normalized_defect_contribution),
                    x.concentration_defect, x.structural_defect)
                    for x in t.complexes if x.name in lu.complexes}) for t in r.results[0].tubes},
            soft_constraints=[c.get() for c in spec.soft_constraints],
            success=r.success)

        conc = Concentrations(
            {T(t): {X(x): TargetConcentration(x.concentration, target=x.target_concentration)
                for x in t.complexes if x.name in lu.complexes} for t in r.results[0].tubes})

        # Build map of undesigned entities to designed entities
        m = Mapping()
        
        for k, v in r.results[0].domains.items():
            m.domains[D(k)] = Domain(v, name=D(k).name, alphabet=D(k).alphabet)

        for s in lu.strands.values():
            m.strands[s] = Strand(''.join(str(m.domains[d]) for d in s.domains), name=s.name)
        for k, x in lu.complexes.items():
            m.complexes[x] = Complex([m.strands[s] for s in x.strands], name=x.name, bonus=x.bonus)
        
        for t in lu.tubes.values():
            m.tubes[t] = Tube({m.strands[s]: v for s, v in t.strand_concentrations().items()},
                SetSpec(0, include=[m.complexes.get(x) or Complex([m.strands[s] for s in x]) for x in t.complexes]),
                name=t.name)

        # lookup for the actual complexes, unnamed -_-
        logq = {x.name: float(x.log_partition_function) for x in r.results[0].complexes}
        pairs = {x.name: x.pair_probabilities.todense()
            for x in r.results[0].complexes if x.pair_probabilities.nnz}

        strands = {s.name: s for s in m.strands.values()}
        complexes = {x.name: Complex([m.strands[s] for s in x.strands], name=x.name, bonus=x.bonus) for x in lu.complexes.values()}

        # Make AnalysisResult with the quantities that have been computed
        a = AnalysisResult()
        for x in r.results[0].complexes:
            x = complexes[x.name]
            pf = Decimal(logq[x.name]).exp() / x.symmetry()
            a.complexes[x] = AnalysisComplexResult(
                pairs=pairs.get(x.name), model=r.model,
                pfunc=pf, free_energy=-1/r.model.beta * float(pf.ln()))

        for t in r.results[0].tubes:
            a.tubes[m.tubes[T(t)]] = AnalysisTubeResult({complexes[x.name]:
                x.concentration for x in t.complexes}, tube=m.tubes[T(t)])

        stats = {k: getattr(r.stats, k) for k in ['num_leaf_evaluations',
            'num_reseeds', 'design_time', 'analysis_time', 'seed']}

        return cls(m, defects, conc, a, stats, r, spec)

    save = Loadable.save
    load = Loadable.load

