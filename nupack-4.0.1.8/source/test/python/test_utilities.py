from nupack import *
import numpy as np

dump = lambda *args: [str(a) for a in args]

################################################################################

def test_strand():
    assert str(Strand('A4', name='a')) == 'AAAA'
    assert str(Strand('AAAA', name='a')) == 'AAAA'
    assert Strand('A4', name='a') == Strand('AAAA', name='a')
    a = Strand('AAAA', name='a')
    assert a == Strand(str(a), name='a')
    a2 = TargetStrand([Domain('AAAA', name='blah')], name='blah2')

def test_domain():
    d = Domain('A4', name='blah')
    e = Domain('A4', name='blah2')
    assert d == d
    assert d != e

################################################################################

def test_pfunc():
    pf, fe = pfunc(['AAAA', 'CCCC'], model=Model())
    assert pf == 0
    assert fe == float('inf')

    pf, fe = pfunc(['AAAA'], model=Model())
    assert pf == 1
    assert float(pf.ln()) == 0
    assert fe == 0

################################################################################

def test_pairs():
    P = pairs(('AAAA', 'TTTT'), model=Model()).to_array()
    assert np.array_equal(P, P.T)
    assert np.all(abs(P.sum(0) - 1) < 1e-8)
    assert np.all(P >= 0)
    assert np.all(P <= 1)

################################################################################

def test_mfe():
    structures = mfe(['A4', 'T4'], model=Model())
    dump(structures)

################################################################################

def test_energy():
    e = energy(['AAAAT'], '(...)', model=Model())
    dump(e)

################################################################################

def test_prob():
    p = structure_probability(['AAAAT'], '(...)', model=Model())
    dump(p)

################################################################################

def test_count():
    n = ensemble_size(['AAAAT'], model=Model())
    dump(n, type(n))

################################################################################

def test_subopt():
    s = subopt(['AAAATTTT'], energy_gap=1.5, model=Model())
    dump(s)

################################################################################

def test_sample():
    s = sample(['AAAATTTT'], num_sample=100, model=Model())
    dump(s)

################################################################################

def test_loop_energy():
    dGloop1 = loop_energy(['AA', 'TT'],
        model=Model(material='RNA', ensemble='stacking')) # stack energy

################################################################################

def test_structure_energy():
    dGstruc1 = structure_energy(['AAAA', 'TTTT'], structure='((((+))))',
        model=Model(material='DNA', celsius = 25, ensemble='stacking'))

################################################################################

def test_design():
    designed_sequences = des('(((+)))', model=Model())
    dump(designed_sequences)
    # --> ['CCC', 'GGG']

# ################################################################################

def test_defect():
    my_defect = defect('(((+)))', ['CCC', 'GGG'], model=Model())
    dump(my_defect)
    # --> ['CCC', 'GGG']

def test_seq_distance():
    assert seq_distance('AA+TT', 'AA+CC') == 2
    assert seq_distance('TT', 'CC') == 2
    assert seq_distance(['TT'], ['CC']) == 2

    assert seq_distance('AA+TT', 'RR+ST') == 1
    assert seq_distance('RR+SS', 'RR+ST') == 1


def test_struc_distance():
    assert struc_distance('.....', '(...)') == 2
    assert struc_distance('.7', '((.3))') == 4
