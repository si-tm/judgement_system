from ..utility import check_instance, check_instances, Loadable, long_output, \
    printable, NamedTuple, from_argument, BaseResult, create_name, mappable, Rules

from ..analysis import TubeResult, Result as AnalysisResult, ComplexResult
from ..core import Structure, Domain, Domain, TargetStrand, TargetComplex, Sequence, TargetTube, SetSpec

from . import results, components
from .core import Options, Specification as Design, Optimization
from .constraints import Match, Pattern, Window, SSM, Library, Diversity, Similarity, Complementarity, EnergyMatch
from .weights import Weights
from .results import Result, ComplexDefect, TubeDefect, TargetDefect, Defects, TargetConcentration, Concentrations

from .components import Timer, TimeInterval, MutationInterval, \
    WriteToFileCheckpoint, StopCondition


################################################################################

tube_design = Design

################################################################################

def complex_design(complexes, model, *, hard_constraints=(),
    soft_constraints=(), defect_weights=None, options=None):
    '''
    Design a set of complexes to adopt specified structures
    - complexes (list[TargetComplex]): complexes to design
    - model (nupack.Model): thermodynamic model to use
    - hard_constraints (list): hard constraints that must be satisfied at all times during the design
    - soft_constraints (list): soft constraints that are added as a weighted penalty to the objective
    - defect_weights  (nupack.design.Weights): specification for tubes, complexes, domains, and strands
    - options (nupack.design.Options): design algorithm options
    A tube is made for each complex, and a Result is returned
    '''
    tubes = [TargetTube({c: 1e-8}, SetSpec(0), name='Tube[%s]' % c.name)
        for c in check_instances(complexes, TargetComplex)]
    return Design(tubes, model, hard_constraints=hard_constraints,
        soft_constraints=soft_constraints, defect_weights=defect_weights, options=options)

################################################################################

def structure_design(structure, strands=None, *, model, options=None):
    '''
    A raw utility meant to take a structure string and return a list of domain sequence strings
    - structure(str or Structure): structure to design
    - strands(List[str]): nucleotide codes of each strand (optional)
    - model(Model): free energy model
    - options(DesignOptions): specific design options
    '''
    structure = Structure(structure)
    if strands is None:
        domains = [Domain('N' * l,  name='Domain[%d]' % i) for i, l in enumerate(structure.lengths())]
    else:
        strands = strands.split('+') if isinstance(strands, str) else strands
        domains = [Domain(s, name='Domain[%d]' % i) for i, s in enumerate(strands)]
    strands = [TargetStrand([d], name='Strand[%d]' % i) for i, d in enumerate(domains)]
    complexes = [TargetComplex(strands, structure=structure, name='complex')]
    design = complex_design(complexes, model=model, options=options)
    result = design.run_one()
    return [str(result.domains[d]) for d in domains]

des = structure_design

################################################################################

def structure_defect(structure, strands, *, model):
    '''
    A raw utility meant to take a structure string and evaluate the multitube defect
    - structure(str or Structure): target structure
    - strands(List[str]): list of strand sequences (fully determined)
    - model(Model): free energy model
    '''
    structure = Structure(structure)
    strands = strands.split('+') if isinstance(strands, str) else strands
    domains = [Domain(str(s),  name='Domain[%d]' % i) for i, s in enumerate(strands)]
    strands = [TargetStrand([d], name='Strand[%d]' % i) for i, d in enumerate(domains)]
    complexes = [TargetComplex(strands, structure=structure, name='complex')]
    design = complex_design(complexes, model=model)
    return design.evaluate().ensemble_defect

defect = structure_defect

################################################################################
