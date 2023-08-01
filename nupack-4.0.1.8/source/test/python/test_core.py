from nupack import *

def test_version():
    import nupack
    assert '@' not in nupack.__version__

################################################################################

def test_strands():
    A = Strand('AGTCTAGGATTCGGCGTGGGTTAA', name='A')
    assert A.nt() == 24
    B = Strand('TTAACCCACGCCGAATCCTAGACTCAAAGTAGTCTAGGATTCGGCGTG', name='B')
    C = Strand('AGTCTAGGATTCGGCGTGGGTTAACACGCCGAATCCTAGACTACTTTG', name='C')

    a1 = Domain('AGTCTAGGATTCGGCGT', name='a1')
    a2 = Domain('GGGTTAA', name='a2')
    D = TargetStrand([a1, a2], name='D') # mostly useful in a design context
    assert D.nt() == 24

    c1 = Complex([A])
    c2 = Complex([A, B, B, C])
    c3 = Complex([A, A])
    c4 = Complex([A, B, C])
    c1a = Complex([A, B])
    c1b = Complex([B, A])
    c1c = Complex([A, B])

    assert c1a == c1b
    assert c1a == c1a
    assert c1a == c1c

def test_parallelism():
    ncpu = core.SharedExecutor().threads()
    config.parallelism = False
    assert config.executor().threads() == 1
    config.parallelism = True
    assert config.executor().threads() == ncpu
    for i in range(5):
        config.threads = i
        assert config.executor().threads() == i if i else ncpu
    config.parallelism = False
    assert config.executor().threads() == 1

################################################################################