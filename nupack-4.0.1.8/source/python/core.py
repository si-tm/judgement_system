from typing import Tuple, List, NamedTuple
import re, json, numpy, collections, concurrent.futures, multiprocessing
from scipy.sparse import csc_matrix
from .utility import str_repr, DUp_to_dpp, check_instances, create_name, complement_name
from .rebind import forward, Callable
from .std import Vector
from . import rotation

##########################################################################

@forward
class JSON:
    '''C++ JSON class'''

    def load(self, string=None):
        '''Load from JSON string'''

    def dump(self, indent=0) -> str:
        '''Return JSON string with given indent'''

    def to_object(self):
        '''Convert self to a python builtin object'''
        return json.loads(self.dump())

    def dump_binary(self) -> numpy.ndarray:
        '''Convert self to a numpy array'''

    def load_binary(self, state):
        '''Load self from a numpy array'''

    @classmethod
    def from_object(cls, obj):
        '''Convert from a python builtin object'''
        out = cls()
        out.load(json.dumps(obj))
        return out

    @classmethod
    def from_path(cls, path):
        '''Load from a file path'''
        out = cls()
        try:
            out.load_file(str(path))
        except RuntimeError as e:
            raise RuntimeError('Failed to load JSON from {} as {}'.format(path, cls.__name__))
        return out

    def __getstate__(self):
        return self.dump_binary()

    def __setstate__(self, state):
        self.move_from(JSON())
        self.load_binary(memoryview(state))

    def __str__(self):
        return self.dump()

##########################################################################

class Pickleable:
    def new_from_json(self, json, _fun_=None):
        '''Load another instance of this class from a nupack.JSON object'''
        return _fun_(json) if isinstance(json, JSON) else _fun_(JSON(json))

    @classmethod
    def from_json(cls, json):
        return cls.__new__(cls).new_from_json(json)

    def to_json(self) -> JSON:
        '''Convert self to a nupack.JSON object'''

    def __getstate__(self):
        return self.to_json()

    def __setstate__(self, state):
        self.move_from(self.new_from_json(state))

##########################################################################

class SimpleType:
    def __str__(self):
        return '{}({})'.format(type(self).__name__,
            ', '.join(k + '=' + repr(getattr(self, k)) for k in self.__annotations__))

    def annotated(self):
        return {k: getattr(self, k) for k in self.__annotations__}

##########################################################################

@forward
class PairList(Pickleable):
    '''
    A container with a list of base pairs `p` such that  `p[i] == j` if
    `i` is paired to `j` and `p[i] == i` if `i` is unpaired
    '''
    values: Vector

    def __init__(self, pairs, _fun_=None):
        if isinstance(pairs, PairList):
            self.copy_from(pairs)
        elif isinstance(pairs, str):
            if 'D' in pairs or 'U' in pairs:
                pairs = DUp_to_dpp(pairs)
            _fun_(self, pairs)
        else:
            _fun_(self, list(pairs)) #  numpy.asarray(pairs, dtype=numpy.uint32)

    def dp(self, nicks=()) -> str:
        '''Return equivalent dot-parens-plus string'''

    def pairlist(self):
        '''Return numpy array of pair list indices'''
        return self.view().copy()

    def view(self):
        return self.values.cast(numpy.ndarray)

    def matrix(self):
        '''Return square structure matrix of 0s and 1s'''
        S = numpy.zeros([len(self)] * 2, dtype=numpy.int32)
        for i, j in enumerate(self):
            S[i, j] = 1
        return S

    structure_matrix = matrix

    def pseudoknots(self) -> List[Tuple[int,int,int,int]]:
        '''Return list of i,j,k,l pseudoknots'''

    def __repr__(self):
        return 'PairList(%r)' % self.view().tolist()

    def __str__(self):
        return str(self.view())

    def __iter__(self):
        return iter(self.pairlist())

    def __getitem__(self, i):
        return self.view()[i]

    def __len__(self):
        '''Number of nucleotides, synonym of nt()'''
        return len(self.view())

    def nt(self):
        '''Number of nucleotides, synonym of len()'''
        return len(self.view())

    def __xor__(self, other) -> int:
        '''Nucleotide distance with other pair list'''

##########################################################################

def struc_distance(structure1, structure2):
    '''Return integer nucleotide distance between two structures'''
    try:
        assert structure1.nicks == structure2.nicks
    except AttributeError:
        pass
    return PairList(structure1) ^ PairList(structure2)

##########################################################################

@forward
class Structure(PairList):
    '''
    Secondary structure including base pairing and nick information
    '''

    def __init__(self, structure, nicks=None, _fun_=None):
        '''Create Structure from either a dot-parens-plus or dpp-rle string representation'''
        if structure is None:
            _fun_(self)
        elif isinstance(structure, Structure):
            self.copy_from(structure)
        elif isinstance(structure, str):
            if 'U' in structure or 'D' in structure:
                structure = DUp_to_dpp(structure)
            _fun_(self, structure)
        else:
            assert nicks is not None
            _fun_(self, structure, list(nicks))

    def lengths(self):
        '''Return lengths of each strand as a numpy array'''
        return numpy.diff((0,) + tuple(self.nicks())).astype(numpy.uint32)

    def nicks(self) -> numpy.ndarray:
        '''A list N of zero-based indices of each base 3′ of a nick between strands'''

    def dp(self, minimum_run=0) -> str:
        '''
        dpp representation of the structure
        minimum_run: min number of characters to run length encode
        minimum_run < 2 means no run length encoding
        '''

    def dotparensplus(self) -> str:
        '''full dpp representation of the structure'''
        return self.dp(0)

    def rle_dotparensplus(self) -> str:
        '''run-length encoded dpp representation of the structure'''
        return self.dp(3)

    def __str__(self):
        return self.dp()

    def __repr__(self):
        return "Structure('%s')" % str(self)

##########################################################################

@forward
class SharedExecutor:
    '''
    Executor for thermodynamic calculation parallelism
    '''
    def __init__(self, threads=0, _fun_=None):
        _fun_(self, int(threads))

    def threads(self) -> int:
        pass

##########################################################################

@forward
class Local:
    '''
    Executor for serial or shared memory parallelism
    '''
    def __init__(self, threads=None, _fun_=None):
        '''Initialize with given number of threads'''
        if threads is None:
            from nupack import config
            _fun_(self, 0 if config.parallelism else 1)
        else:
            _fun_(self, int(threads))

    def n_workers(self) -> int:
        '''Return number of workers'''

##########################################################################

# def set_sequence_type(rna: int) -> int:
#     '''
#     Set whether sequence should print as RNA; return previous value
#     0: print as RNA (weak control)
#     1: print as RNA (strong control)
#     2: print as DNA (weak control)
#     3: print as DNA (strong control)
#     '''

# class MaterialContext:
#     def __init__(self, rna):
#         self.rna = 1 if rna else 3
#         self.backup = 0

#     def __enter__(self):
#         self.backup = set_sequence_type(self.rna)

#     def __exit__(self, v, t, tb):
#         set_sequence_type(self.backup)

##########################################################################

@forward
class Base(Pickleable):
    '''
    Class for a single base: ACGT or a wildcard of those bases
    U is represented the same as T
    '''
    def __init__(self, integer):
        pass

    def __pos__(self) -> int:
        '''Return integer of this base'''

    def __str__(self):
        return str(+self)

    def __repr__(self):
        return 'Base(%d)' % +self

##########################################################################

@forward
class Wildcard(Pickleable):
    '''
    Class for a single base: ACGT or a wildcard of those bases
    U is represented the same as T
    '''
    def __init__(self, integer_or_base):
        pass

    def __pos__(self) -> int:
        '''Return integer mask of this wildcard'''

    def indices(self) -> List[int]:
        '''Indices of possible bases'''

    def first(self) -> Base:
        '''First possible base in this wildcard'''

    def is_determined(self) -> bool:
        '''Return whether exactly one base is possible'''
    
    def includes(self, base_or_wildcard) -> bool:
        '''Return whether another base specification is a subset of this one'''

    def __xor__(self, other):
        '''Logical mask operator'''
    
    def __and__(self, other):
        '''Logical mask operator'''
    
    def __or__(self, other):
        '''Logical mask operator'''

    def __str__(self):
        return str(+self)

    def __repr__(self):
        return 'Base(%d)' % +self

##########################################################################

def are_compatible(bases1, bases2) -> bool:
    '''
    Return whether two bases, wildcards, sequences, or domains are compatible
    Compatible means that there is a possible sequence simultaneously satisfying both sides of constraints
    '''

##########################################################################

@forward
class Alphabet(Pickleable):
    '''
    Class representing an alphabet of different bases, possibly for multiple materials
    '''

    def __init__(self, string='', name=None):
        '''Deduce an alphabet from a user given sequence string and optional alphabet name'''
    
    def complement(self, base_or_wildcard):
        pass 
    
    def sequence(self, sequence_string):
        pass 
    
    def domain(self, sequence_string):
        pass 
    
    def letter(self, base_or_wildcard) -> str:
        pass 
    
    def reverse_complement(self, sequence):
        pass 

    def length(self) -> int:
        '''Number of bases'''

##########################################################################

def standardize_alphabet(complexes, alphabet):
    '''Return a list of SequenceList from the named complexes given'''

##########################################################################

class _SequenceTraits:

    def __getitem__(self, index, _fun_=None) -> Base:
        '''Return a given base'''
        if abs(index) >= len(self):
            raise IndexError(index)
        return _fun_(self, index % len(self))

    def __contains__(self, base) -> bool:
        '''Test if a given base is contained'''

    def __iter__(self):
        return iter(map(Base, str(self)))

    def __len__(self) -> int:
        '''Number of nucleotides'''

    def nt(self) -> int:
        '''Number of nucleotides'''

##########################################################################

@forward
class Sequence(Pickleable, _SequenceTraits):
    '''
    Class representing a sequence of fixed bases
    '''
    def __init__(self, sequence_or_string, material=None, alphabet=None, _fun_=None):
        if isinstance(sequence_or_string, str):
            if alphabet is None:
                alphabet = Alphabet(sequence_or_string, name=material)
            _fun_(self, sequence_or_string, alphabet)
        else:
            _fun_(self, sequence_or_string)

@forward
class RawDomain(Pickleable, _SequenceTraits):
    '''
    Class representing a sequence of variable bases which may be denoted
    by IUPAC wildcardds
    '''
    def __init__(self, sequence_or_string, material=None, alphabet=None, _fun_=None):
        if isinstance(sequence_or_string, str):
            if alphabet is None:
                alphabet = Alphabet(sequence_or_string, name=material)
            _fun_(self, sequence_or_string, alphabet)
        else:
            _fun_(self, sequence_or_string)

##########################################################################

class _NamedTraits:
    def __invert__(self):
        '''Watson Crick reverse complement'''
        return self.reverse_complement(False)

    def to_string(self, minimum_run=0) -> str:
        '''
        String version of the sequence
        If minimum_run < 2, return raw sequence
        Otherwise, run length encode each run of at least minimum_run
        '''

    def __str__(self):
        return self.to_string(0)

##########################################################################

@forward
class Strand(Sequence, _NamedTraits):
    '''
    Class for an RNA or DNA sequence, possibly containing wildcards
    '''
    name: str
    def __init__(self, string, name, material=None, alphabet=None, _fun_=None):
        '''Initialize from sequence string, name, and (optional) alphabet'''
        if alphabet is None:
            alphabet = Alphabet(string, name=material)
        _fun_(self, Sequence(string, material=material, alphabet=alphabet), create_name(name), alphabet)

    def __repr__(self):
        return '<Strand %s>' % self.name

    def reverse_complement(self, wobble=False):
        '''Watson Crick or wobble reverse complement'''

##########################################################################

@forward
class Domain(RawDomain, _NamedTraits):
    '''
    Class for a domain of RNA or DNA, uniquely named and usually containing wildcards
    '''
    name: str

    def __init__(self, string, name, material=None, alphabet=None, _fun_=None):
        '''Initialize from sequence string, name, and (optional) alphabet'''
        if not isinstance(string, str) and material is None and alphabet is None:
            raise TypeError('No way to deduce alphabet for Domain')
        if alphabet is None:
            alphabet = Alphabet(string, name=material)
        _fun_(self, RawDomain(string, alphabet=alphabet), create_name(name), alphabet)

    def reverse_complement(self, wobble=False):
        '''Return the reverse complement of self'''

    def __add__(self, other):
        '''Return a DomainList with self and other as contents'''

    def __repr__(self):
        return '<Domain %s>' % self.name

    def substitute(self, rules):
        return rules.get(self, self)

##########################################################################

def seq_distance(complex1, complex2, alphabet=None, _fun_=None):
    '''Hamming distance between 2 sequences or complexes'''
    try:
        a = SequenceList(complex1, alphabet=alphabet)
    except Exception:
        a = RawDomainList(complex1, alphabet=alphabet)

    try:
        b = SequenceList(complex2, alphabet=alphabet)
    except Exception:
        b = RawDomainList(complex2, alphabet=alphabet)

    return _fun_(a, b).cast(int)

##########################################################################

@forward
class DomainList(Vector):
    '''
    Utility class for working with a list of domains
    '''

    def __init__(self, list):
        pass

    def nt(self) -> int:
        '''Number of nucleotides'''

    def __add__(self, other):
        '''Return concatenation of self with domains of other'''

    def __repr__(self):
        return 'DomainList(%r)' % list(self)

    def __str__(self) -> str:
        '''String representation'''
        return '[%s]' % ', '.join(d.name for d in self)

    def __contains__(self, domain) -> bool:
        '''Test if a domain is contained'''

    def reverse_complement(self, wobble=False):
        '''Return the reverse complement of self'''

    def __invert__(self):
        '''Watson Crick reverse complement'''
        return self.reverse_complement(False)

##########################################################################

@forward
class TargetStrand(Strand):
    '''
    Class for a designable strand of RNA or DNA, uniquely named and containing Domains
    '''
    name: str
    domains: Tuple[Domain, ...]

    def __init__(self, domains, *, name, _fun_=None):
        '''Initialize from list of domains'''
        domains = check_instances(domains, Domain)
        _fun_(self, domains, create_name(name))

    def __iter__(self):
        return iter(self.domains)

    def __len__(self):
        '''Number of domains'''
        return len(self.domains)

    def ndomains(self):
        '''Number of domains'''
        return len(self.domains)

    def __getitem__(self, index):
        return self.domains[index]

    def __str__(self):
        return self.name

    def __repr__(self):
        return '<TargetStrand %s>' % self.name

    def reverse_complement(self, wobble=False):
        '''Return the reverse complement of self'''

    def substitute(self, rules):
        return rules.get(self) or TargetStrand([d.substitute(rules) for d in self.domains], name=self.name)

##########################################################################

@forward
class SequenceListBase(Pickleable):
    def __getitem__(self, index, _fun_=None) -> Strand:
        if abs(index) >= len(self):
            raise IndexError(index)
        return _fun_(self, index % len(self))

    def __len__(self) -> int:
        pass

    def lowest_rotation(self):
        '''Return equivalent complex to self, but rotated to canonical ordering'''

    def nicks(self) -> List[int]:
        '''List of strand break positions'''

    def symmetry(self) -> int:
        '''Rotational symmetry factor'''

    def nt(self) -> int:
        '''Number of nucleotides'''

##########################################################################

@forward
class RawDomainList(SequenceListBase):
    def __init__(self, string_or_strings, alphabet=None, _fun_=None):
        v = string_or_strings
        v = re.split('[+, ]+', v) if isinstance(v, str) else v
        _fun_(self, [RawDomain(s, alphabet=alphabet) for s in v])

##########################################################################

@forward
class SequenceList(SequenceListBase):
    def __init__(self, string_or_strings, alphabet=None, _fun_=None):
        v = string_or_strings
        v = re.split('[+, ]+', v) if isinstance(v, str) else v
        _fun_(self, [Sequence(s, alphabet=alphabet) for s in v])

    # @property
    # def strands(self):
    #     return tuple(self)

    # def nstrands(self):
    #     '''Number of strands'''
    #     return len(self)

    # def __iter__(self):
    #     return (self[i] for i in range(len(self)))

    # def __contains__(self, strand) -> bool:
    #     '''Test whether a given strand is contained'''

    # def __str__(self):
    #     return '+'.join(map(str, self))

    # def __repr__(self):
    #     return 'Complex([%s])' % ', '.join(map(repr, map(str, self)))

##########################################################################

@forward
class Complex(Pickleable):
    '''
    Class for an ordered complex of RNA or DNA, not containing wildcards
    '''
    strands: Tuple[Strand, ...]
    bonus: float
    name: str

    def __init__(self, strands, name='', bonus=0.0, _fun_=None):
        '''
        Initialize from a list of strands
        - strands: list of Sequence or Strand, or SequenceList or TargetComplex
        '''
        _fun_(self, check_instances(strands, Strand), create_name(name), float(bonus))

    def lowest_rotation(self):
        '''Return equivalent complex to self, but rotated to canonical ordering'''
        
    #     assert not isinstance(strands, (Sequence, Strand)), 'strands should be a list of strands. Try using SequenceList([strand])'
    #     if isinstance(strands, TargetComplex):
    #         _fun_(self, list(strands))
    #     else:
    #         _fun_(self, strands)

    def nt(self) -> int:
        '''Return number of nucleotides'''

    def nstrands(self) -> int:
        '''Return number of strands'''

    def __iter__(self):
        return iter(self.strands)

    def __len__(self):
        return len(self.strands)

    def nicks(self) -> List[int]:
        '''List of strand break positions'''

    def __repr__(self):
        return '<Complex %s>' % self.name
        
    # def __xor__(self, other) -> int:
    #     '''Hamming distance in sequences'''


    # def __getstate__(self):
    #     return str(self)

    # def __setstate__(self, state):
    #     self.copy_from(SequenceList(state))

##########################################################################

@forward
class TargetComplex(Pickleable):
    '''
    Class for an ordered complex of RNA or DNA, containing designable strands
    and possibly an associated target structure
    '''
    name: str
    structure: Structure
    bonus: float
    strands: Tuple[TargetStrand, ...]

    def __init__(self, strands, structure=None, *, bonus=0, name=None, _fun_=None):
        '''
        strands: list of TargetStrand, or another TargetComplex
        '''
        strands = check_instances(strands, TargetStrand)
        name = create_name(name) if name else ''
        if structure is None or isinstance(structure, str):
            structure = Structure(structure)
        if not isinstance(structure, Structure):
            raise TypeError('Expected Structure or str: got %r' % structure)
        _fun_(self, strands, structure, name, float(bonus))

    def lowest_rotation(self):
        '''Return equivalent complex to self, but rotated to canonical ordering'''

    def nstrands(self):
        '''Number of strands'''
        return len(self)

    def __len__(self) -> int:
        return len(self.strands)

    def __getitem__(self, index):
        return self.strands[index]

    def __iter__(self):
        return iter(self.strands)

    # def __contains__(self, strand) -> bool:
    #     '''Test whehther a given TargetStrand is contained'''

    def nt(self) -> int:
        '''Number of nucleotides'''

    def __str__(self):
        return '+'.join(map(str, self.strands))

    def __repr__(self):
        return '<TargetComplex %s>' % self.name

##########################################################################

@forward
def rotational_symmetry(indices) -> int:
    '''Rotational symmetry of a list of indices'''

@forward
def compute_necklaces(callback: Callable[[List[int]], None], *, size, n_elements) -> int:
    '''compute necklaces and pass them to callback function'''


def complexes_from_max_size(strands, max_size):
    complexes = set()

    def add_perm(l):
        names = tuple(rotation.lowest_rotation([strands[i] for i in l]))
        complexes.add(names)

    for i in range(1, max_size+1):
        compute_necklaces(add_perm, n_elements=len(strands), size=i)

    return complexes

##########################################################################

class SetSpec(NamedTuple):
    '''
    Specification for a set of complexes that will be generated given a set of strands
    '''
    max_size: int = 1
    include: tuple = ()
    exclude: tuple = ()

    def generate(self, strands, cls=tuple):
        '''Generate a list of complexes from a given set of strands and complex type'''
        complexes = set(map(cls, self.include))
        for c in complexes_from_max_size(strands, self.max_size):
            complexes.add(cls(c))
        complexes.difference_update(map(cls, self.exclude))
        return sorted(complexes)

    @classmethod
    def get(cls, obj):
        if obj is None:
            return cls()
        if isinstance(obj, cls):
            return obj
        return cls(0, obj)

##########################################################################

@forward
class ComplexSet(Pickleable):
    '''
    Class representing a set of complexes with shared strands
    '''
    strands: List[Strand]
    complexes: List[Complex]

    def __init__(self, strands, complexes=SetSpec(), _fun_=None):
        strands = check_instances(strands, Strand)
        spec = SetSpec.get(complexes)
        complexes = spec.generate(strands, Complex)
        _fun_(self, strands, complexes)

    def __iter__(self):
        return iter(self.complexes)

    # @classmethod
    # def union(cls, complex_sets):
    #     '''Set union over an iterable of instances of ComplexSet'''
    #     return cls(complexes=[c for t in complex_sets for c in t.complexes])

    # def __add__(self, other):
    #     '''Set union with another instance (see .union())'''
    #     return self.union([self, other])

    def __str__(self):
        return 'ComplexSet(%s, include=%s)' % (self.strands, self.complexes)

    def __repr__(self):
        return 'ComplexSet(%r, include=%r)' % (self.strands, self.complexes)

################################################################################

@forward
class Tube(ComplexSet):
    '''
    Class representing a set of complexes in a tube with optional strand concentrations
    '''
    name: str
    concentrations: numpy.ndarray

    def __init__(self, strands, complexes=SetSpec(), *, name, _fun_=None):
        name = create_name(name)
        strands, conc = zip(*dict(strands).items())
        conc = numpy.array(conc, dtype=numpy.float64)
        _fun_(self, ComplexSet(strands, complexes), conc, name)

    # @classmethod
    # def union(cls, tubes):
    #     '''
    #     Set union over an iterable of instances of Tubee
    #     Concentrations of the same strands are added together
    #     '''
    #     d = {}
    #     for t in tubes:
    #         for s, c in zip(t.strands, t.concentrations):
    #             d[s] = d.get(s, 0) + c
    #     return cls(d, complexes=[c for t in tubes for c in t.complexes],
    #         name='+'.join(t.name for t in tubes))

    def __str__(self):
        return 'Tube({%s}, name=%r)' % (', '.join('{}: {}'.format(s.name, c)
            for s, c in zip(self.strands, self.concentrations)), self.name)

    def __repr__(self):
        return '<Tube %s>' % self.name
    #     return 'Tube(%r, complexes=%r, name=%r)' % (dict(zip(self.strands,
    #         self.concentrations)), list(self.complexes), self.name)

##########################################################################

# @mappable
@forward
class TargetTube(Pickleable):
    name: str
    complexes: Tuple[TargetComplex, ...]
    concentrations: numpy.ndarray
    n_on_targets: int

    def __init__(self, on_targets, off_targets=SetSpec(), *, name, _fun_=None):
        '''
        Initialize tube from:
        - on_targets: a dict of TargetComplex to concentration (float)
        - off_targets: a dict containing optional keys:
            - include: an additional iterable of TargetComplex to include as offtargets
            - exclude: an additional iterable of TargetComplex to exclude
            - max_size (default 1): if >0, all non-on-target complexes of up to this size will be included as offtargets
        - name: name of the tube (str)
        '''
        on_targets = dict(on_targets)
        name = create_name(name)

        for t, conc in on_targets.items():
            if not t.structure:
                raise ValueError('Target complex %r was not given a structure' % t)
            assert conc >= 0, 'Concentration must be positive'

        concentrations = numpy.array(list(on_targets.values()), dtype=numpy.float64)

        spec = SetSpec.get(off_targets)
        strands = tuple(set(s for v in (on_targets, spec.include, spec.exclude) for x in v for s in x))
        off_targets = spec.generate(strands, TargetComplex)

        complexes = list(on_targets) + list(off_targets)
        _fun_(self, complexes, concentrations, name)

    def strand_concentrations(self):
        '''dict of strands to their concentrations'''
        out = {}
        for k, v in zip(self.on_targets, self.concentrations):
            for s in k.strands:
                out[s] = out.get(s, 0) + v
        return out

    @property
    def on_targets(self):
        return self.complexes[:self.n_on_targets]

    @property
    def off_targets(self):
        return self.complexes[self.n_on_targets:]

    def __str__(self):
        return 'TargetTube({%s}, name=%r)' % (', '.join('%s: %.2e' % (k.name, v)
            for k, v in zip(self.on_targets, self.concentrations)), self.name)

    def __repr__(self):
        return '<TargetTube %s>' % self.name

    # @property
    # def complexes(self):
    #     '''All complexes in the tube (read-only)'''
    #     return self.off_targets.union(self.on_targets)

    # def replace(self, rules):
    #     self.off_targets = rules.map(self.off_targets)
    #     self.strands = rules.map(self.strands)
    #     self.on_targets = {rules(k): v for k, v in self.on_targets.items()}

    # def __lt__(self, other):
    #     if isinstance(other, TargetTube):
    #         return self.name < other.name
    #     return NotImplemented

    # def __eq__(self, other):
    #     if not isinstance(other, TargetTube):
    #         return NotImplemented
    #     return (self.name, self.strands, self.on_targets, self.off_targets) == \
    #         (other.name, other.strands, other.on_targets, other.off_targets) and \
    #         numpy.array_equal(self.concentrations, other.concentrations)

    # def __hash__(self):
    #     return hash(self.name)

##########################################################################

@forward
class Fenwick:
    def __len__(self):
        pass

##########################################################################

@forward
class MemoryLimit:
    length: int
    capacity: int

##########################################################################

@forward
class LRU:
    pass

################################################################################

@forward
class Sparsity(Pickleable):
    '''
    Specification of a sparsity level for pair probability calculation
    - row_size: maximum number of non-diagonal elements per row to include
    - threshold: the minimum pair probability to be considered
    '''
    threshold: float
    row_size: int

    def __init__(self, *, threshold=0, row_size=0):
        pass

################################################################################

@forward
class PairMatrix(Pickleable):
    '''
    Raw minimal representation of a sparse pair probability matrix
    - diagonal: value of each diagonal element
    - values: value of each offdiagonal entry (i, j) s.t. i < j
    - rows: row i of each offdiagonal entry
    - cols: column j of each offdiagonal entry
    '''
    diagonal: numpy.ndarray
    values: numpy.ndarray
    rows: numpy.ndarray
    cols: numpy.ndarray

    def __init__(self, full, *, fraction, threshold, _fun_=None):
        '''Initialize from full matrix and desired sparsity between 0 and 1'''
        if fraction < 0 or fraction > 1:
            raise ValueError('Sparsity fraction should be in [0:1] (is {})'.format(f))
        if threshold < 0 or threshold > 1:
            raise ValueError('Sparsity threshold should be in [0:1] (is {})'.format(t))

        if hasattr(full, 'tocsc'):
            assert full.shape[0] == full.shape[1], full.shape
            assert fraction == 0, 'fraction should be 0 for sparse input'
            assert threshold == 0, 'threshold should be 0 for sparse input'
            full = full.tocsc().astype(numpy.float64)
            _fun_(self, full)
        else:
            full = numpy.asfortranarray(full, dtype=numpy.float64)
            assert full.ndim == 2, full.shape
            assert full.shape[0] == full.shape[1], full.shape
            _fun_(self, full, Sparsity(threshold=threshold, row_size=round(full.shape[0] * fraction)))

    def defect(self, structure):
        '''Calculate structure defect with respect to a Structure or PairList'''
        P = self.to_array()
        return len(P) - (structure.matrix() * P).sum()

    def to_array(self) -> numpy.ndarray:
        '''Convert to a numpy array'''

    def to_sparse(self, cls=csc_matrix):
        '''Return equivalent sparse matrix'''
        n = len(self.diagonal)
        rows = numpy.concatenate([numpy.arange(n), self.rows, self.cols])
        cols = numpy.concatenate([numpy.arange(n), self.cols, self.rows])
        vals = numpy.concatenate([self.diagonal, self.values, self.values])
        return cls((vals, (rows, cols)), shape=(n, n))

    def __str__(self):
        with numpy.printoptions(formatter={'float': '{:.4f}'.format}):
            return str(self.to_array())

################################################################################

class _Config:
    @property
    def parallelism(self):
        return self._parallelism

    @parallelism.setter
    def parallelism(self, value):
        assert isinstance(value, bool)
        self._parallelism = value
        self.threads = 0 if value else 1

    @property
    def threads(self):
        return self._threads

    @threads.setter
    def threads(self, value):
        from .core import SharedExecutor
        self._threads = int(value)
        self._parallelism = (self._threads != 1)
        self._executor = SharedExecutor(self._threads)

    def executor(self):
        if self._executor is None:
            self.threads = self.threads
        return self._executor

    @classmethod
    def augment(cls, other):
        type(other).executor = cls.executor
        type(other).threads = cls.threads
        type(other).parallelism = cls.parallelism
        type(other)._cpu_count = multiprocessing.cpu_count()
        assert other._cpu_count > 0
        other.thread_pool = concurrent.futures.ThreadPoolExecutor(other._cpu_count)