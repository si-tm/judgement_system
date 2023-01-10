from functools import reduce
import numpy as np
import subprocess
import os
import re
import scipy.stats
import math


def convert_domain_list(pepperDomainList):
    res = []
    for p in pepperDomainList:
        res.append((p.name,p.sequence))
        pp= p.getComplementary()
        res.append((str(pp),p.getCompSequence()))
    return res

def domainsdict_addcompl(dom):
    pepperd = {a: PepperDomain(a, 5, sequence=dom[a]) for a,d in dom.items()}
    return {**dom, **{f"{a}*":p.getCompSequence() for a,p in pepperd.items()} }


class Reaction:
    def __init__(self, idReac, reactants, products, rate, r_type = None):
        self.idReac = idReac
        self.reactants = reactants
        self.products = products
        self.rate = rate
        self.r_type = r_type

#Parse a reaction from PepperCorn
def parseReaction(trimmedLine):
    reactionRate = float(re.search(r'(?<==)\s+[-\w.+]+', trimmedLine).group(0).strip())
    idReac = re.search(r'(?<=])\s+[^\n]+', trimmedLine).group(0).strip()
    species = idReac.split('->')
    reactants = [a.strip() for a in species[0].split('+')]
    products = [a.strip() for a in species[1].split('+')]
    r_type = re.search(r'(?<=\[)(\w|-)+', trimmedLine).group(0).strip()
    return Reaction(idReac, reactants, products, reactionRate, r_type)

def basePairing(dna):
    if dna in {'T', 't'}:
        return 'A'
    elif dna in {'A', 'a'}:
        return 'T'
    elif dna in {'G', 'g'}:
        return 'C'
    elif dna in {'C', 'c'}:
        return 'G'
    return 'X'

#DNA domains
class PepperDomain(object):

    def __init__(self, name, length, compl=False, sequence=""):
        self.name = name
        self.length = length #Used for structure stability calculation in PepperCorn, independent from actual sequence length
        self.compl = compl
        self.sequence = sequence

    def getComplementary(self):
        return PepperDomain(self.name,self.length, compl=not self.compl, sequence = self.sequence)

    def getCompSequence(self):
        return ''.join([basePairing(s) for s in self.sequence[::-1]])

    def getPepperDescr(self):
        return "length "+str(self.name)+" = "+str(self.length)

    def __str__(self):
        return str(self.name)+('*' if self.compl else '')

#Standard domains from the Kakenhi project
domaina = PepperDomain('a',5, sequence="TCCCT")
domainb = PepperDomain('b',5, sequence="GAGTC")

#DNA strands, made of domains or str
class PepperDNAStrand(object):

    def __init__(self, name, domains, concentration = 0.0):
        self.name = name
        self.domains = domains #must be iterable
        self.concentration = concentration

    def getPepperDescr(self):
        return str(self.name)+" = "+reduce(lambda x,y: str(x)+" "+str(y),self.domains)+(" @ initial "+str(self.concentration)+" M" if self.concentration > 0.0 else "")

    #Takes a domain to sequence mapping and return the DNA sequence of the strand
    def getFullSequence(self):
        return ''.join([a.sequence if not a.compl else a.getCompSequence() for a in self.domains])

    def __str__(self):
        return "PepperDNAStrand("+str(self.name)+", "+str([str(a) for a in self.domains])+")"

#DNA systems, defines domains and strands
class PepperDNASystem(object):

    def __init__(self, domains, strands):
        self.domains = domains #must be iterable; should not contain complementary
        self.strands = strands #must be iterable

    def getPepperDescr(self):
        s = reduce(lambda x,y: str(x)+"\n"+str(y),[a.getPepperDescr() for a in self.domains])+"\n"
        return s+ reduce(lambda x,y: str(x)+"\n"+str(y),[a.getPepperDescr() for a in self.strands])

    def getFullSequences(self):
        return [s.getFullSequence() for s in self.strands]

    def __str__(self):
        return "PepperDNASystem("+str([str(a) for a in self.domains])+","+str([str(a) for a in self.strands])+")"

def getAllPossibleDomains(baseDomains):
    return [a for a in baseDomains]+[a.getComplementary() for a in baseDomains]

#For genome evaluation, get a PepperDNAStrand from a value
def generateStrandFromID(name, value, fullDomains, length = 4):
    doms = [ ]
    #print(value)
    tmpVal,conc = value
    nDomains = len(fullDomains)
    for i in range(length):
        doms.append(fullDomains[int(tmpVal % nDomains)])
        tmpVal = tmpVal / nDomains
    doms.reverse()
    return PepperDNAStrand(name,doms,concentration=conc)

#For genome evaluation, get a system from a genome
def generateDNASystem(values,domains = [domaina,domainb], length = 4):
    #print(values)
    doms = getAllPossibleDomains(domains)
    strands = [generateStrandFromID('s'+str(i),v,doms,length = length) for i,v in enumerate(values)]
    return PepperDNASystem(domains,strands)


#Function to submit a system
def submitSystem(dnaSystem, evalName = "test.pil", outputFile = "result.pil", logFile = "pepper.log", maxComplexSize=200, maxComplexCount=200, maxReactionCount=1000, BFS=True):
    with open(evalName,'w') as inpFile:
        inpFile.write(dnaSystem.getPepperDescr())
    executable = ['peppercorn', '-o',outputFile,evalName,"--max-complex-size=%i" % maxComplexSize, "--max-complex-count=%i" % maxComplexCount, "--max-reaction-count=%i" % maxReactionCount]
    if BFS:
        executable.append('--bfs-ish')
    if logFile is not None:
        with open(logFile, 'w') as f:
            #subprocess.run(['peppercorn', '-o',outputFile,evalName,'--max-complex-size=200', '--max-complex-count=20000', '--max-reaction-count=100000'], stdout=f, stderr=subprocess.STDOUT)
            subprocess.run(executable, stdout=f, stderr=subprocess.STDOUT)
    else:
        subprocess.run(executable, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

def getEmptyResult():
    res = {'length' : 0, 'largestStruct': 0, 'smallestStruct' : 0, 'medianStruct': 0, 'meanStruct': 0, 'mostCommonSize': 0, 'structSizeStddev': 0,
            'totalConnections': 0, 'medianConnections': 0, 'connectionsVariance' : 0, 'reactionCount': 0, 'log10ReactionCount': 0, 'reactionRatio': 0, 'invReactionRatio': 0, 'structSizeQuartileCoeffDispersion': 0, 'structSizeMAD': 0, 'structSizeNormalisedStddev': 0, 'structSizeCoeffVar': 0, 'entropyComplexSize': 0, 'nbLoops': 0, "skewComplexSize": 0, "nbModesComplexSize": 0, "entropyReactionTypes": 0, "bindReactionFreq": 0, "branchReactionFreq": 0, "nbBindReactions": 0, "nbBranchReactions": 0, "log10nbBindReactions": 0, "log10nbBranchReactions": 0}
    original = list(res.keys())
    #for k in original:
    #    res['delta'+k] = 0
    return res

#Could be replaced by counting every items, then looking at the biggest one
#That would be advantageous if we switch to extra features like that in the future
def getMaxMostCommonSize(li):
    countsSoFar = {}
    mostCommonID = -1
    for i in li:
        if i in countsSoFar :
            countsSoFar[i] += 1
        else :
            countsSoFar[i] = 1
        if mostCommonID == -1 or countsSoFar[i] > countsSoFar[mostCommonID]:
            mostCommonID = i
    return mostCommonID

# DEPRECATED
#Parse a reaction
#def parseReaction(trimmedLine):
#    reactionRate = float(re.search(r'(?<==)\s+[-\w.+]+', trimmedLine).group(0).strip())
#    species = re.search(r'(?<=])\s+[^\n]+', trimmedLine).group(0).strip().split('->')
#    return {'reactants' : [a.strip() for a in species[0].split('+')], 'products': [a.strip() for a in species[1].split('+')], 'rate' : reactionRate}


_reactions_types = ["condensed", "open", "bind11", "bind21", "branch-3way", "branch-4way"]


#Function to analyse the results from evaluation
def getPropsFile(pilFile, logFile, domainSize = 5, maxComplexSize = 200):
    res = getEmptyResult()
    try:
        with open(pilFile) as f:
            allText = f.read()
        m = re.search('(?<=Resting complexes \n)[^#]+',allText)
        allStabs = m.group(0).strip().split('\n')
        allSizes = [s.count('+')+1 for s in allStabs]
        allConnections = [s.count('(') for s in allStabs]
        allSizes.sort()
        allConnections.sort()
        largest = max(allSizes)
        totalConnections = reduce(lambda x,y : x + y, allConnections)
        m = re.search('(?<=reactions).+',allText,flags=re.DOTALL)
        #print(m.group(0))
        allReactions = m.group(0).strip().split('\n')

        reactions = [parseReaction(line) for line in allReactions]
        reactions_distrib = {key: len([r for r in reactions if r.r_type == key]) for key in _reactions_types}
        reactions_tot = sum(reactions_distrib.values())
        reactions_freq = [v / reactions_tot for v in reactions_distrib.values()]
        res['entropyReactionTypes'] = scipy.stats.entropy(reactions_freq)
        res['nbBindReactions'] = reactions_distrib['bind11'] + reactions_distrib['bind21']
        res['nbBranchReactions'] = reactions_distrib['branch-3way'] + reactions_distrib['branch-4way']
        if res['nbBindReactions'] > 0:
            res['log10nbBindReactions'] = math.log10(res['nbBindReactions'])
        else:
            res['log10nbBindReactions'] = 0.
        if res['nbBranchReactions'] > 0:
            res['log10nbBranchReactions'] = math.log10(res['nbBranchReactions'])
        else:
            res['log10nbBranchReactions'] = 0.
        res['bindReactionFreq'] = reactions_freq[2] + reactions_freq[3]
        res['branchReactionFreq'] = reactions_freq[4] + reactions_freq[5]

        #print(allReactions)
        res['length'] =  len(allStabs)
        res['largestStruct'] =  largest
        res['smallestStruct'] = min(allSizes)
        res['totalConnections'] = domainSize*totalConnections
        res['mostCommonSize'] = getMaxMostCommonSize(allSizes)
        res['medianStruct'] = np.median(allSizes)
        res['meanStruct'] = np.mean(allSizes)
        res['medianConnections'] = np.median(allConnections)
        res['structSizeStddev'] = np.std(allSizes)
        res['reactionCount'] = len(allReactions)
        res['log10ReactionCount'] = math.log10(len(allReactions)) if len(allReactions) else 0.
        res['reactionRatio'] = res['length'] / res['reactionCount'] #res['reactionCount'] / res['length']
        res['invReactionRatio'] = res['reactionCount'] / res['length']
        q1, q3 = np.percentile(allSizes, [25, 75], interpolation="linear")
        res['structSizeQuartileCoeffDispersion'] = (q3 - q1) / (q3 + q1)
        #print("TMP:", res['structSizeQuartileCoeffDispersion'], q1, q3)
        res['structSizeMAD'] = np.median([abs(x - np.median(allSizes)) for x in allSizes]) # Mean Absolute deviation
        #print("TMPMAD:", res['structSizeMAD'])
        res['structSizeNormalisedStddev'] = np.std(allSizes) / max(allSizes)
        res['structSizeCoeffVar'] = np.std(allSizes) / np.mean(allSizes)

        allSizes = np.array(allSizes)
        possibleComplexSizes = list(range(1, maxComplexSize+1))
        nbComplexesOfSize = [np.sum(allSizes == x) for x in possibleComplexSizes]
        freqComplexesOfSize = np.array(nbComplexesOfSize) / float(len(allSizes))
        res['entropyComplexSize'] = scipy.stats.entropy(freqComplexesOfSize)
        res['skewComplexSize'] = scipy.stats.skew(freqComplexesOfSize)
        res['nbModesComplexSize'] = len(scipy.stats.mode(freqComplexesOfSize)[0])

        if logFile is not None:
            with open(logFile) as f:
                allLog = f.read()
            res['nbLoops'] = allLog.count("WARNING: Double stem count in Loop() Object.")

    except:
        print("WARNING: Peppercorn didn't work out")
        raise
    return res


class ChemicalSystem:
    def __init__(self, defConc):
        self.species = {}
        self.reactions = {}
        self.speciesSizes = {}
        self.speciesConnections = {}
        self.speciesReactions = {}
        self.maximumOrderRates = {}
        self.defaultConc = defConc


def setMaximumOrderReactionRates(system):
    for r in system.reactions.values():
        se = set(r.reactants)
        for s in se:
            count = r.reactants.count(s)
            maxOrder = system.maximumOrderRates[s] if s in system.maximumOrderRates else ()
            if not maxOrder or (maxOrder[0] < count or (maxOrder[0] == count and maxOrder[1] < r.rate)):
                system.maximumOrderRates[s] = (count,r.rate)


#Parse a pil file generated by PepperCorn
def getSystemFromPil(pilfile, defaultConc = 1e-7):
    system = ChemicalSystem(defaultConc)
    with open(pilfile,'r') as pil:
        for line in pil:
            l = line.strip()
            if len(l) == 0 or l[:6] == "length" or l[0] == '#':
                pass
            elif l[:8] == "reaction":
                r = parseReaction(l)
                setOfReactants = set(r.reactants)
                for s in setOfReactants:
                    if s in system.speciesReactions:
                        system.speciesReactions[s].append(r)
                    else:
                        system.speciesReactions[s] = [r]
                system.reactions[r.idReac] = r
            else:
                # species
                speciesData = l.split('=')
                name = speciesData[0].strip()
                size = speciesData[1].count('+') + 1
                connections = speciesData[1].count('(')
                conc = float(re.search(r'(?<=initial)\s+[-\w.+]+', l).group(0).strip()) if '@' in l else 0.0
                system.species[name] = int(conc/defaultConc)
                system.speciesSizes[name] = size
                system.speciesConnections[name] = connections
    setMaximumOrderReactionRates(system)
    return system


#Top-level function, length is the length in domain of the DNA strands
def evaluateGenotype(genotype, domains = [domaina,domainb], length = 4, configName = "test.pil", outputName = "result.pil", logName = "pepper.log", maxComplexSize=200, maxComplexCount=200, maxReactionCount=1000, BFS=True):
    #print(repr(genotype))
    #print("DEBUG1: ", configName, outputName, logName)
    system = generateDNASystem([(i,x) for i, x in enumerate(genotype) if x > 0.0], domains = domains, length = length)
    if not system.strands:
        return getEmptyResult()
    #Just in case peppercorn bugs, write a dummy output file
    submitSystem(system, evalName = configName, outputFile = outputName, logFile = logName, maxComplexSize = maxComplexSize, maxComplexCount = maxComplexCount, maxReactionCount = maxReactionCount, BFS = BFS)
    return getPropsFile(outputName, logName, maxComplexSize=maxComplexSize)

def getSequenceListFromGenotype(genotype, domains = [domaina,domainb], length = 4):
    
    system = generateDNASystem([(i,x) for i, x in enumerate(genotype) if x > 0.0], domains = domains, length = length)
    return system.getFullSequences()

# MODELINE	"{{{1
# vim:expandtab:softtabstop=4:shiftwidth=4:fileencoding=utf-8
# vim:foldmethod=marker
