import yaml

from . import Model, Design, TargetTube, TargetComplex, EnergyMatch, \
    Match, Complementarity, Similarity, Library, Window, Pattern, Diversity, SSM

################################################################################

class BlockList(list):
    pass

def list_yaml(dumper, data):
    return dumper.represent_sequence(u'tag:yaml.org,2002:seq', data, flow_style=False)

yaml.add_representer(BlockList, list_yaml)

################################################################################

class BlockDict(dict):
    pass

def dict_yaml(dumper, data):
    return dumper.represent_mapping( u'tag:yaml.org,2002:map', data, flow_style=False )

yaml.add_representer(BlockDict, dict_yaml)

################################################################################

def model_yaml(dumper, model):
    o = {
        'temperature': model.temperature,
        'sodium': model.conditions.na_molarity,
        'magnesium': model.conditions.mg_molarity,
        'material': model.parameters.info.file.path,
        'ensemble': str(model.ensemble)
    }
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', o, flow_style=False)

yaml.add_representer(Model, model_yaml)

################################################################################

def target_tube_yaml(dumper, tube):
    f = lambda v: [x.name if x.is_named() else [s.name for s in x.strands]
                  for x in map(TargetComplex, v)]
    o = {
        'strands': {s.name: float(v) for s, v in zip(tube.strands, tube.concentrations)},
        'name': tube.name,
        'max_size': tube.spec.max_size,
        'include': f(tube.spec.include),
        'exclude': f(tube.spec.exclude)
    }
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', o, flow_style=False)

yaml.add_representer(TargetTube, target_tube_yaml)

################################################################################

def match_yaml(dumper, c):
    o = {'regions': c.regions}
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', o, flow_style=False)

yaml.add_representer(Match, match_yaml)

################################################################################

def complementarity_yaml(dumper, c):
    o = {'regions': c.regions, 'wobble_mutations': c.wobble_mutations}
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', o, flow_style=False)

yaml.add_representer(Complementarity, complementarity_yaml)

################################################################################

def similarity_yaml(dumper, c):
    o = {
        'domains': c.domains,
        'reference': c.reference,
        'limits': c.limits,
        'weight': c.weight,
    }
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', o, flow_style=False)

yaml.add_representer(Similarity, similarity_yaml)

################################################################################

def library_yaml(dumper, c):
    o = {
        'domains': c.domains,
        'catalog': c.catalog,
    }
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', o, flow_style=False)

yaml.add_representer(Library, library_yaml)

################################################################################

def window_yaml(dumper, c):
    o = {
        'domains': c.domains,
        'sources': c.sources,
    }
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', o, flow_style=False)

yaml.add_representer(Window, window_yaml)

################################################################################

def pattern_yaml(dumper, c):
    o = {
        'pattern': pattern,
        'scope': scope,
        'weight': weight,
    }
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', o, flow_style=False)

yaml.add_representer(Pattern, pattern_yaml)

################################################################################

def diversity_yaml(dumper, c):
    o = {}
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', o, flow_style=False)

yaml.add_representer(Diversity, diversity_yaml)

################################################################################

def ssm_yaml(dumper, c):
    o = {}
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', o, flow_style=False)

yaml.add_representer(SSM, ssm_yaml)

################################################################################

def energymatch_yaml(dumper, c):
    o = {}
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', o, flow_style=False)

yaml.add_representer(EnergyMatch, energymatch_yaml)

################################################################################

def design_yaml(dumper, design):
    '''Return a YAML string of the design specification'''
    o = {}
    o['domains'] = BlockDict({d.name: d.to_string(4) for d in design.domains})
    # o['domains'] = [{'name': d.name, 'sequence': str(d)} for d in design.domains]
    o['strands'] = {s.name: [d.name for d in s.domains] for s in design.strands}

    o['complexes'] = {c.name: {'strands': [s.name for s in c.strands],
        'structure': c.structure.dp(4) if len(c.structure) else None}
        for c in design.complexes if len(c.structure) or c.is_named()}

    o['tubes'] = BlockList(design.tubes)

    o['model'] = design.model

    if design.hard_constraints:
        o['hard_constraints'] = BlockList(design.hard_constraints),

    if design.soft_constraints:
        o['soft_constraints'] = BlockList(design.soft_constraints),

    o['options'] = design.options.to_dict()

    x = lambda c: [c.name if c.is_named() else [s.name for s in c.strands]]
    o['defect_weights'] = [dict(domain=d.name, strand=s.name, complex=x(c), tube=t.name, weight=w)
        for _, (d, s, c, t, w) in design.defect_weights.iterrows() if w != 1]

    if design.objective_weight != 1:
        o['objective_weight'] = design.objective_weight

    return dumper.represent_mapping(u'tag:yaml.org,2002:map', o, flow_style=False)

yaml.add_representer(Design, design_yaml)

################################################################################
