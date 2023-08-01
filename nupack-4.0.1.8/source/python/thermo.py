from .rebind import Tuple, List, Dict, Callable, forward
from .core import Structure, LRU, PairList, Complex, PairMatrix, Local, Sparsity
from .model import Model, Ensemble
from .utility import match, check_instances

import numpy, decimal, typing

################################################################################

class Job:
    '''Variant class holding a complex and a specified job input (e.g. PFJob)'''
    def __init__(self, complex, job):
        pass

################################################################################

class PFJob:
    '''Low-level class representing a partition function job'''
    pass

################################################################################

class MFEJob:
    '''Low-level class representing an MFE energy job'''
    pass

################################################################################

class PairsJob:
    '''Low-level class representing a pair probabilities job'''
    sparsity: Sparsity

    def __init__(self, sparsity):
        pass

################################################################################

class CostsJob:
    '''Low-level class representing a MFE costs job'''

################################################################################

class SampleJob:
    '''
    Low-level class representing a Boltzmann sampling job.
    Note that `seed` is the real seed used for the random number generation.
    That is, seed=0 has no special meaning.
    '''
    number: int
    seed: int

    def __init__(self, number=1, seed=0, _fun_=None):
        seed = seed or numpy.random.randint(0, 1 << 63, dtype=numpy.int64)
        assert seed >= 0, 'Random number generation seed should be non-negative'
        _fun_(self, number, seed)

################################################################################

class SuboptJob:
    '''Low-level class representing a suboptimal structures job'''
    gap: float
    max_number: int

    def __init__(self, gap=0.0, max_number=100000):
        pass

################################################################################

class Result:
    '''
    Variant class holding optional specific results for a given complex
    - pfunc:  optional PFResult
    - mfe:    optional MFEResult
    - pairs:  optional PairsResult
    - sample: optional SampleResult
    - subopt: optional SuboptResult
    '''

class PFResult:
    '''Result for partition function calculation'''
    logq: float
    raw_logq: float

class MFEResult:
    '''Result for MFE energy calculation'''
    energy: float
    raw_energy: float

class PairsResult:
    '''Result for pair probability calculation'''

class CostsResult:
    '''Result for MFE costs calculation'''

class SampleResult:
    '''Result for Boltzmann sample calculation'''

class SuboptResult:
    '''Result for suboptimal structure calculation'''

################################################################################

class Future:
    '''Asynchronous future for a thermo job submission'''

    def get(self):
        '''Wait for the result and retrieve it'''

################################################################################

class ComputeOptions:
    '''Cache and execution options for analysis jobs'''
    max_bytes: int

    def __init__(self, *, max_bytes, _fun_=None):
        '''Initialize from number of threads and max_bytes'''
        from nupack import config
        if max_bytes < 0:
            max_bytes = int(config.cache * 1e9)
        assert max_bytes >= 0
        _fun_(self, config.executor(), max_bytes)

################################################################################

@forward
def submit(jobs, model, options, _fun_=None):
    '''Submit thermo jobs and return a Future for their result'''
    jobs = check_instances(jobs, Job)
    assert isinstance(model, Model)
    assert isinstance(options, ComputeOptions)
    return _fun_(jobs, model, options, gil=False)

################################################################################

@forward
class PairsEnergy:
    '''Class holding pairlist and associated energy for a given state'''
    structure: PairList
    energy: float
    stack_energy: float

################################################################################

class StructureEnergy(typing.NamedTuple):
    '''Class holding structure and associated energy for a given state'''
    structure: Structure
    energy: float
    stack_energy: float

    def __repr__(self):
        return 'StructureEnergy(%r, energy=%s, stack_energy=%s)' % (self.structure, self.energy, self.stack_energy)

    def __str__(self):
        return 'StructureEnergy(%r, energy=%.2f, stack_energy=%.2f)' % (str(self.structure), self.energy, self.stack_energy)

################################################################################
