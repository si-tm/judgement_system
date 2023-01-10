import pickle
from glob import glob
import sys
from qdpy.phenotype import Individual

class Custom_unpickler(pickle._Unpickler):

    current_module = {"Individual": "qdpy.phenotype", "Fitness": "qdpy.phenotype", "IndividualLike": "qdpy.phenotype", }
    
    def find_class(self, module, name):
        if name in Custom_unpickler.current_module:
            module = Custom_unpickler.current_module[name] #backward compatibility
        sys.audit('pickle.find_class', module, name)
        if self.proto < 3 and self.fix_imports:
            if (module, name) in _compat_pickle.NAME_MAPPING:
                module, name = _compat_pickle.NAME_MAPPING[(module, name)]
            elif module in _compat_pickle.IMPORT_MAPPING:
                module = _compat_pickle.IMPORT_MAPPING[module]
        __import__(module, level=0)
        if self.proto >= 4:
            return pickle._getattribute(sys.modules[module], name)[0]

        else:
            return getattr(sys.modules[module], name)

def get_res(filename):
    res = None
    with open(filename,"rb") as f:
        res = Custom_unpickler(f,fix_imports=True, encoding="ASCII", errors="strict",
          buffers=None).load()
    return res

p = '../../input/p/seqA-GA100000-0.80/final_20200904131038.p'
results = get_res(p)

#print(results["container"].best_index)

#(1, 28, 37)

# ind = results["container"].solutions[results["container"].best_index][0]
# #[Individual([0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0])]
# _domains = {'a': 13, 'b': 13}
# dnaDomains = [submitPepperCorn.PepperDomain(a,_domains[a]) for a in _domains.keys()]
# length = 2

# #outputFile = "~/"
# outputFile = "/Users/takepy/takeoxdna/kakenhievolvedna/scripts/output"

#  # Evaluate individual
# configFile = outputFile + "conf.pil"
# pilFile = outputFile + "output.pil"
# #logFile = outputFile + "log.pil"
# with open (configFile,"w") as f:
#     f.write("test")
# scores = submitPepperCorn.evaluateGenotype(ind, domains=dnaDomains, 
# length=length, maxComplexSize=30, 
# maxComplexCount=400, maxReactionCount=2000, 
# configName=configFile, outputName=pilFile, 
# logName=logFile)


#def evaluateGenotype(genotype, domains = [domaina,domainb], length = 4, configName = "test.pil", outputName = "result.pil", logName = "pepper.log", maxComplexSize=200, maxComplexCount=200, maxReactionCount=1000, BFS=True):

#result = submitPepperCorn.evaluateGenotype(***, domains = [T,T], length = 4, configName = "test.pil", outputName = "result.pil", logName = "pepper.log", maxComplexSize = 200, maxComplexCount=200, maxReactionCount=1000, BFS=True)

#2.pilファイルを得る