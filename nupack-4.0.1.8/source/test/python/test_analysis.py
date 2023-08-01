from nupack import *
import numpy as np, pickle

################################################################################

def pickles(a):
    aname = getattr(a, 'name', None)
    s = pickle.dumps(a)
    b = pickle.loads(s)
    bname = getattr(b, 'name', None)
    assert a == b
    assert aname == bname

def test_pickle():
    A = Strand('AGTCTAGGATTCGGCGTGGGTTAA', name='A')
    pickles(A)
    a1 = Domain('AGTCTAGGATTCGGCGT', name='a1')
    pickles(a1)
    D = TargetStrand([a1], name='D')
    pickles(D)
    c1 = SequenceList([A])
    pickles(c1)
    c1 = Complex([A], name='c1')
    pickles(c1)
    tc = TargetComplex([D], structure='.17', name='A')
    pickles(tc)

def test_strands():
    #DNA tests
    a_d = Model(material = 'dna04').alphabet()
    A = Strand('AGTCTAGGATTCGGCGTGGGTTAA', name='A', alphabet=a_d)
    B = Strand('TTAACCCACGCCGAATCCTAGACTCAAAGTAGTCTAGGATTCGGCGTG', name='B', alphabet=a_d)
    C = Strand('AGTCTAGGATTCGGCGTGGGTTAACACGCCGAATCCTAGACTACTTTG', name='C', alphabet=a_d)
    A_d = Strand('dTTT', name='A_d', alphabet=a_d)
    try:
        print(Strand('rTTT', name='r_in_d', alphabet=a_d))
        assert False
    except RuntimeError:
        pass

    #RNA tests
    a_r = Model().alphabet()
    A_r = Strand('AAA', name='A_r', alphabet=a_r)
    B_r = Strand('TTT', name='B_r', alphabet=a_r)
    C_r = Strand('rTTT', name='C_r', alphabet=a_r)
    try:
        print(Strand('dTTT', name='d_in_r', alphabet=a_r))
        assert False
    except RuntimeError:
        pass

    a1 = Domain('AGTCTAGGATTCGGCGT', name='a1')
    a2 = Domain('GGGTTAA', name='a2')
    D = TargetStrand([a1, a2], name='D') # mostly useful in a design context

    c1 = SequenceList([A])
    c2 = SequenceList([A, B, B, C])
    c3 = SequenceList([A, A])
    c4 = SequenceList([A, B, C])
    c1a = SequenceList([A, B])
    c1b = SequenceList([B, A])
    c1c = SequenceList([A, B])

    assert c1a == c1b
    assert c1a == c1a
    assert c1a == c1c

def test_sample():
    a = Strand('C10G10', name='a')
    b = Strand('C10G10', name='b')
    cs = ComplexSet([a, b], complexes=SetSpec(max_size=2))
    model = Model()
    res = complex_analysis(cs, compute=['sample'], model=model, options={'num_sample': 5})
    for k, v in res.complexes.items():
        assert len(v.sample) == 5

################################################################################

def test_named():
    A1 = Strand('AGTCTAGGATTCGGCGTGGGTTAA', name='A1')
    A2 = Strand('AGTCTAGGATTCGGCGTGGGTTAA', name='A2')

    assert A1 != A2 # --> False
    assert A1 == A1 # --> True

    c2 = Complex([A1, A2])
    c1 = Complex([A1, A1])

    t1 = Tube({A1: 1e-6, A2: 1e-8}, name='')
    t2 = Tube({A1: 1e-6, A2: 1e-9}, complexes=SetSpec(3, include=[c2], exclude=[c1]), name='')

################################################################################

def test_tube_analysis():
    #DNA
    m = Model(material = 'dna04')
    a = m.alphabet()
    A = Strand('CTGATCGAT', name='A', alphabet = a)
    B = Strand('GATCGTAGTC', name='B', alphabet = a)

    t1 = Tube({A: 1e-8, B: 1e-9}, complexes=SetSpec(3), name='0')
    t2 = Tube({A: 1e-10, B: 1e-9}, complexes=SetSpec(2), name='1')

    result = tube_analysis(tubes=[t1, t2],
        compute=['pairs', 'mfe', 'sample'], model=m,
        options={'num_sample': 100})

    t1 = Tube({A: 1e-8, B: 1e-8}, complexes=SetSpec(2), name='t1')
    result = complex_analysis(t1 + t2, compute=['pairs', 'mfe'], model=m)
    result[Complex([A])] # --> ComplexResult

    t1_result = complex_concentrations(t1, result, concentrations={A: 1e-8, B: 1e-9}) # use manually specified concentrations if desired
    t2_result = complex_concentrations(t2, result) # use concentration from t2

    _ = str(t1_result) # result concentrations

    #RNA
    m = Model()
    a = m.alphabet()
    A = Strand('CTGATCGAT', name='A', alphabet = a)
    B = Strand('GATCGTAGTC', name='B', alphabet = a)
    
    t1 = Tube({A: 1e-8, B: 1e-9}, complexes=SetSpec(3), name='0')
    t2 = Tube({A: 1e-10, B: 1e-9}, complexes=SetSpec(2), name='1')

    result = tube_analysis(tubes=[t1, t2],
        compute=['pairs', 'mfe', 'sample'], model=m,
        options={'num_sample': 100})

    t1 = Tube({A: 1e-8, B: 1e-8}, complexes=SetSpec(2), name='t1')
    result = complex_analysis(t1 + t2, compute=['pairs', 'mfe'], model=m)
    result[Complex([A])] # --> ComplexResult

    t1_result = complex_concentrations(t1, result, concentrations={A: 1e-8, B: 1e-9}) # use manually specified concentrations if desired
    t2_result = complex_concentrations(t2, result) # use concentration from t2

    _ = str(t1_result) # result concentrations

################################################################################

def test_tube_analysis_4():
    #DNA
    m = Model(material = 'dna04')
    alpha = m.alphabet()
    a = Strand('CAGTCGATC', name='a', alphabet = alpha)
    b = Strand('ATCGACGTA', name='b', alphabet = alpha)
    c = Complex([a, b])

    t1 = Tube({a: 1e-6, b: 1e-9}, complexes=SetSpec(include=[c]), name='t1')
    t2 = Tube({a: 1e-8, b: 1e-9}, complexes=SetSpec(include=[c]), name='t2')

    result = tube_analysis([t1, t2], compute=['pairs', 'mfe', 'sample', 'subopt'],
        options={'num_sample': 2, 'energy_gap': 0.5}, model=m)
    s = str(result)
    s = result._repr_html_()

    result[t1] # --> TubeResult

    result[t1].complex_concentrations # --> [1.5e-10]
    result[t1].ensemble_pair_fractions # --> [[1.0, 0.0], [0.0, 1.0]]

    result[c]
    # pfunc for complex c
    result[c].pfunc
    # mfe for complex c
    list(result[c].mfe)
    # ppairs matrix for complex c
    result[c].pairs

    #RNA
    m = Model()
    alpha = m.alphabet()
    a = Strand('CAGTCGATC', name='a', alphabet = alpha)
    b = Strand('ATCGACGTA', name='b', alphabet = alpha)
    c = Complex([a, b])

    t1 = Tube({a: 1e-6, b: 1e-9}, complexes=SetSpec(include=[c]), name='t1')
    t2 = Tube({a: 1e-8, b: 1e-9}, complexes=SetSpec(include=[c]), name='t2')

    result = tube_analysis([t1, t2], compute=['pairs', 'mfe', 'sample', 'subopt'],
        options={'num_sample': 2, 'energy_gap': 0.5}, model=m)
    s = str(result)
    s = result._repr_html_()

    result[t1] # --> TubeResult

    result[t1].complex_concentrations # --> [1.5e-10]
    result[t1].ensemble_pair_fractions # --> [[1.0, 0.0], [0.0, 1.0]]

    result[c]
    # pfunc for complex c
    result[c].pfunc
    # mfe for complex c
    list(result[c].mfe)
    # ppairs matrix for complex c
    result[c].pairs

################################################################################

def test_partition_function():
    pf, fe = pfunc('A100', model=Model(ensemble='some-nupack3'))

    assert pf == 1
    assert fe == 0

    #DNA
    m = Model(material = 'dna04')
    pf, fe = pfunc('dA100', model=m)
    assert pf == 1
    assert fe == 0
    #RNA
    m = Model()
    pf, fe = pfunc('rA100', model=m)
    assert pf == 1
    assert fe == 0


################################################################################

def test_mfe():
    #DNA
    m = Model(material = 'dna04')
    mf = mfe('dA100', model=m)
    assert mf[0].energy == 0
    assert mf[0].structure == Structure('.' * 100)
    #RNA
    m = Model()
    mf = mfe('rA100', model=m)
    assert mf[0].energy == 0
    assert mf[0].structure == Structure('.' * 100)

################################################################################

def test_ensemble_size():
    assert ensemble_size('GAAAC', model=Model()) == 2
    assert ensemble_size('GAAA', model=Model()) == 1

################################################################################

def test_impossible_complex():
    strands = 'A+A'
    model = Model(ensemble='some-nupack3')
    pf, fe = pfunc(strands, model=model)
    samples = sample(strands, model=model, num_sample=10)
    mfes = mfe(strands, model=model)
    subopts = subopt(strands, energy_gap=2, model=model)
    ps = pairs(strands, model=model)

    assert pf == 0
    assert fe == np.inf
    assert np.all(ps.to_array() == np.eye(2))
    assert len(mfes) == 0
    assert len(subopts) == 0
    assert len(samples) == 10

################################################################################

def test_overflow():
    model = Model(ensemble='some-nupack3', material='rna95-nupack3')
    pf, _ = pfunc('C300G300', model=model)
    assert abs(float(pf.ln()) - 1393.4600830078125) < 0.01

    P = pairs('C300G300', model=model)
    assert abs(P.to_array().sum(0) - 1).max() < 0.01
    assert abs(P.to_sparse().sum(0) - 1).max() < 0.01

################################################################################

def test_unpaired():
    t = Tube({Strand('AAAA', name='A'): 1e-8}, name='t')
    res = tube_analysis([t], model=Model()).tubes[t]
    assert res.fraction_bases_unpaired is None
    res = tube_analysis([t], model=Model(), compute=['pairs']).tubes[t]
    assert res.fraction_bases_unpaired == 1.0

