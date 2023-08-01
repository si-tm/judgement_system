from typing import Tuple, List, NamedTuple
import json

# nupack includes
from .. import rotation  # import rotation.lowest_rotation
from .. import constants  # import compute_necklaces
from ..model import Model

from ..core import TargetTube, TargetComplex, Domain, Pickleable, JSON
from ..utility import check_instance, Loadable, check_instances
from .results import Result, RawResult
from .weights import Weights

################################################################################

try:
    from .gui import GraphicOptimization as Optimization
except ImportError:
    from .trials import Optimization

################################################################################

class Env:
    def __init__(self, _fun_=None):
        from .. import config
        _fun_(self, config.executor())

################################################################################

class Options:
    '''Algorithm options for design'''

    seed: int = 0
    f_stop: float = 0.02
    f_passive: float = 0.01
    H_split: int = -1 # defaults to 2 for RNA, 3 for DNA
    N_split: int = 12
    f_split: float = 0.99
    f_stringent: float = 0.99
    dG_clamp: float = -20
    M_bad: int = 300
    M_reseed: int = 50
    M_reopt: int = 3
    f_redecomp: float = 0.03
    f_refocus: float = 0.03
    f_sparse: float = 1e-05
    max_time: float = 0
    slowdown: float = 0
    log: str = ''
    decomposition_log: str = ''
    thermo_log: str = ''
    time_analysis: bool = False
    wobble_mutations: bool = False

    def __init__(self, options=None, _fun_=None, **kws):
        _fun_(self, options) if isinstance(options, Options) else _fun_(self)

        if isinstance(options, dict):
            kws.update(options)

        for k, v in kws.items():
            if k not in self.__annotations__:
                raise KeyError('Invalid keyword %r' % k)
            try:
                setattr(self, k, v)
            except TypeError as e:
                raise TypeError('Incorrect type for keyword %r' % k) from e

################################################################################

class Lookup(NamedTuple):
    domains: dict
    strands: dict
    complexes: dict
    tubes: dict

################################################################################

class Specification(Loadable, Pickleable):
    '''
    A design for a set of tubes
    - tubes (list[TargetTube]): tubes to design
    - model (nupack.Model): thermodynamic model to use
    - hard_constraints (list): hard constraints that must be satisfied at all times during the design
    - soft_constraints (list): soft constraints that are added as a weighted penalty to the objective
    - defect_weights  (nupack.design.Weights): specification for tubes, complexes, domains, and strands
    - options (nupack.design.Options): design algorithm options
    '''
    tubes: Tuple[TargetTube, ...]
    model: Model
    options: Options
    objective_weight: float

    def __init__(self, tubes, model, options=None, hard_constraints=(), 
                 soft_constraints=(), defect_weights=None, objective_weight=1.0, _fun_=None):

        tubes = check_instances(tubes, TargetTube)
        options = Options(options)
        defect_weights = Weights(tubes) if defect_weights is None else defect_weights
        model = check_instance(model, Model)
        
        _fun_(self, tubes, model, options, hard_constraints, soft_constraints, 
              defect_weights.factors() if hasattr(defect_weights, 'factors') else defect_weights, objective_weight)

    def replace(self, domains, _fun_=None):
        domains = tuple(check_instances(domains, Domain))
        _fun_(self, domains + tuple(~d for d in domains))

    def apply(self, domains=None, tubes=None, model=None, hard_constraints=None,
        soft_constraints=None, defect_weights=None, objective_weight=None, options=None):
        '''Apply new attributes to a design, returning a new one'''
        new = Specification(
            tubes=self.tubes if tubes is None else tubes,
            model=self.model if model is None else model,
            soft_constraints=self.soft_constraints if soft_constraints is None else soft_constraints,
            hard_constraints=self.hard_constraints if hard_constraints is None else hard_constraints,
            defect_weights=self.defect_weights if defect_weights is None else defect_weights,
            objective_weight=self.objective_weight if objective_weight is None else objective_weight,
            options=self.options if options is None else options,
        )
        if domains is not None:
            new.replace(domains)
        return new

    def run_one(self, restart=None, *, checkpoint_condition=None, checkpoint_handler=None, max_bytes=0, seed_offset=0, _fun_=None):
        """create the designer from the fully specified design and run and
        return results

        Arguments:
            restart: a Result object that can be used to initialize the
                sequence state of the design.
            checkpoint_condition: A 2-argument callable that takes the current design's
                statistics and timer and outputs True if a checkpoint should be emitted
                and False otherwise, e.g. if time has elapsed a certain amount
            checkpoint_handler: A single argument callable that receives a
                Result object and decides how to process it as a checkpoint,
                e.g. printing to a file

        Returns: Result: the result of a finished design
        """
        if restart is not None:
            restart = restart.raw
        if self.options.seed:
            self = self.copy()
            self.options.seed += seed_offset
        try:
            raw = _fun_(self, Env(), checkpoint_condition, checkpoint_handler, restart, gil=False)
        except RuntimeError as e:
            if 'Space::clone' in str(e) or 'Failed constraint space' in str(e):
                raise RuntimeError('Failed to clone the constraint space. This likely means '
                    'your design specification is impossible due to 1) base-pairing or sequence '
                    'length constraints inherent in the design or 2) user-specified constraints. '
                    'The best approach is currently to check your specification or simplify your '
                    'design until it works, allowing you to pinpoint the problem')
            raise
        return Result.build(self, raw)

    def run(self, trials, restart=None, *, checkpoint_condition=None, checkpoint_handler=None):
        r, c, h = ([None] * trials if x is None else list(x) for x in
            (restart, checkpoint_condition, checkpoint_handler))
        assert len(r) == trials, 'Incorrect number of restarts'
        assert len(c) == trials, 'Incorrect number of checkpoint conditions'
        assert len(h) == trials, 'Incorrect number of checkpoint handlers'

        from .. import config
        self.options.cache_bytes_of_RAM = int(config.cache * 1e9 / trials)

        if trials == 1:
            return [self.run_one(restart=r[0], checkpoint_condition=c[0], checkpoint_handler=h[0])]

        tasks = [config.thread_pool.submit(self.run_one, R, checkpoint_condition=C,
            checkpoint_handler=H, seed_offset=i) for i, (R, C, H) in enumerate(zip(r, c, h))]
        return [t.result() for t in tasks]

    def launch(self, trials, *, checkpoint=None, restart=None, interval=600):
        '''
        Launch a series of trials of the design optimization in a non-blocking manner
        - trials: the number of trials
        - checkpoint: the directory to save results (defaults to nowhere)
        - restart: the directory to load starting sequences (defaults to nowhere)
        - interval: how often to request a root evaluation to save results to a file (in seconds)
        '''
        from .. import config
        self.options.cache_bytes_of_RAM = int(config.cache * 1e9 / trials)

        out = Optimization(self, trials=trials, checkpoint=checkpoint, interval=interval)
        out.start(restart)
        return out

    def redo_complement(self, d):
        '''Ugly method to apply custom wobble pair option to deduce the intended complement sequence'''
        assert d.name[-1] == '*'
        try:
            orig = next(x for x in self.domains if x.name == d.name[:-1])
            out = orig.reverse_complement(self.options.wobble_mutations)
            return out
        except StopIteration:
            return d

    def evaluate(self, _fun_=None):
        """
        Compute the multistate test tube ensemble defect for a set of fixed
        sequences or throw and exception if the input sequences have
        variable nucleotides. Returns a Result object.
        """
        return Result.build(self, _fun_(self, Env(), gil=False))

    def load_result(self, path):
        '''Load a Result object from a JSON file path'''
        return Result.build(self, RawResult(json=JSON.from_path(path)))

    def lookup(self):
        domains = {}
        strands = {}
        complexes = {}
        tubes = {}
        for t in self.tubes:
            tubes[t.name] = t
            for c in t.complexes:
                complexes[c.name] = c
                for s in c.strands:
                    strands[s.name] = s
                    for d in s.domains:
                        domains[d.name] = d
        return Lookup(domains, strands, complexes, tubes)

    def to_yaml(self) -> str:
        '''Return YAML string of self'''

    @staticmethod
    def from_yaml(string):
        '''Load a design from YAML string'''



