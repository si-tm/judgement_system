from ..utility import check_instances, check_instance
from ..core import Domain, RawDomain, Sequence, TargetComplex

################################################################################

class Compare:
    def __init__(self, left, right, _fun_=None):
        '''
        Initialize from matching list of domains:
        - left: first list of domains
        - right: second list of domains of equal total length
        '''
        _fun_(self, check_instances(left, Domain), check_instances(right, Domain))

################################################################################

class Match(Compare):
    '''Add constraint making left and right lists of domains identical

    Args:
        left (list[Domain]): list of domains to constrain
        right (list[Domain]): list of domains to constrain (of equal total length)
    '''

################################################################################

class Complementarity:
    '''
    Add constraint making left and right lists of domains complementary
    '''
    def __init__(self, left, right, wobble_mutations=False, _fun_=None):
        '''
        Initialize from matching list of domains:
        - left (list[Domain]): list of domains to constrain
        - right (list[Domain]): list of domains to constrain (of equal total length)
        - wobble_mutations (bool): whether to allow wobble complements
        '''
        _fun_(self, check_instances(left, Domain), check_instances(right, Domain), bool(wobble_mutations))

################################################################################

class Pairing(Compare):
    '''Add constraint making left and right lists of domains complementary

    Args:
        left (list[Domain]): list of domains to constrain
        right (list[Domain]): list of domains to constrain (of equal total length)
    '''

################################################################################

class Similarity:
    '''Make a design element match some fraction of a reference sequence

    Args:
        domains (list[Domain]): list of domains to constrain
        reference (str): a degenerate base sequence
        limits (tuple(float, float)): min and max fraction that must be matched
        weight (float): weight to apply, if using as soft constraint
    '''

    soft_kind = 'similarity'
    soft_name = 'Soft constraints: similarity'

    def __init__(self, domains, reference, *, limits, weight=1, _fun_=None):
        limits = (0, 1) if limits is None else tuple(map(float, limits))
        assert len(limits) == 2
        _fun_(self, check_instances(domains, Domain), RawDomain(reference), limits, float(weight))

    def replace(self, rules):
        self.domains = rules.map(self.domains)

################################################################################

class Library:
    '''Constrain the concatenation of domains to be drawn from some
            alternative concatenation of library sequences

    Args:

        domains list[Domain]: domains to be concatenated in the constraint
        catalog (list[list[Sequence]]): the list of libraries whose sequences will be applied to the constraint
    '''
    def __init__(self, domains, catalog, _fun_=None):
        domains = getattr(domains, 'domains', domains)
        _fun_(self, check_instances(domains, Domain), [[RawDomain(s) for s in c] for c in catalog])

    def replace(self, rules):
        self.domains = rules.map(self.domains)

################################################################################

class Window:
    '''Constrain the concatenation of domains to be drawn as a window
        from one of the sources

    Args:
        domains (list[Domain]): concatenated domains to constrain
        sources (list[Sequence]): list of possible source sequences

    Raises:
        ValueError: if one of the sources listed is shorter than the
            concatenation of domains and hence has no windows.
    '''
    def __init__(self, domains, sources, _fun_=None):
        _fun_(self, check_instances(domains, Domain), [RawDomain(s) for s in sources])

    def replace(self, rules):
        self.domains = rules.map(self.domains)

################################################################################

class Pattern:
    '''Constrain given domains to exclude specified patterns

    Args:
        patterns (list[RawDomain]): list of pattern sequences
        scope (list[Domain]): optional region to apply the constraint (global if not specified)
        weight (float): weight to use if using as soft constraint
    '''

    soft_kind = 'pattern'
    soft_name = 'Soft constraints: pattern'

    def __init__(self, patterns, *, weight=1, scope=None, _fun_=None):
        _fun_(self, [] if scope is None else check_instances(scope, Domain),
            [RawDomain(p) for p in patterns], float(weight))

################################################################################

class Diversity:
    '''Constrain the minimum diversity (minimum_nucleotide_types) per
    windows of word_length within each domain/strand in names

    Args:
        word (int): word length to consider
        diversity (int): number of different nucleotide types that must be present in each word
        scope (list[Domain]): optional region to apply the constraint (global if not specified)
    '''

    def __init__(self, word, types, *, scope=None, _fun_=None):
        _fun_(self, [] if scope is None else check_instances(scope, Domain),
              int(word), int(types))

################################################################################

class SSM:
    '''Degree of sequence symmetry violation for a set of complexes

    Args:
        complexes (list[TargetComplex] or None): list of complexes to constrain (global if None)
        word (int): SSM word size
        weight (float): weight to apply to the SSM objective
    '''

    soft_kind = 'sequence_symmetry'
    soft_name = 'Soft constraints: sequence symmetry'

    def __init__(self, *, scope=None, word=4, weight=1.0, _fun_=None):
        _fun_(self, () if scope is None else check_instances(scope, TargetComplex), int(word), float(weight))

    def replace(self, rules):
        self.scope = self.scope and rules.map(self.scope)

################################################################################

class EnergyMatch:
    '''Constrain difference in duplex energy between sets of domains

    Args:
        domains (list[Domain]): list of domains to constrain
        energy_ref (float): constrain duplex energies toward this reference (=mean if not given)
        weight (float): weight to apply to the objective
    '''

    soft_kind = 'energy difference'
    soft_name = 'Soft constraints: energy difference'

    def __init__(self, domains, *, energy_ref=None, weight=1.0):
        _fun_(self, check_instances(domains, Domain), 
            None if energy_ref is None else float(energy_ref), float(weight))

    def replace(self, rules):
        self.domains = rules.map(self.domains)
