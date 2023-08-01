#!/usr/bin/env python
# coding: utf-8

from qdpy import algorithms, containers, plots
import math
import numpy as np
from keras.models import load_model

import tensorflow as tf
print("TensorFlow version:", tf.__version__)
#tf.get_logger().setLevel('ERROR')

import nupack
from nupack import *
from nupack import Model as nupackModel
from tensorflow.keras import Model
from scipy.stats import linregress
import numpy as np

from qdpy.phenotype import Individual, GenIndividuals
from typing import Iterable, Any
import numpy as np
from numpy import random

from qdpy import algorithms, containers, plots, tools
from qdpy.containers import Container
from qdpy.algorithms.evolution import Evolution
from qdpy.base import DomainLike, registry
# from qdpy import registry
import math
from functools import partial
import functools
import csv

class L3Individual(Individual):
    strands: Iterable[tuple[Strand,float]]
    complexes: ComplexSet
    evaluation: Any
    # input/structure_seq/input_seq_L3.csv
    # [0, 1, 1, 0 ... 1] : a a*, a b, ... b* b*
    indexes: [bool] 
    structure : [str] 
    domains : {}

    def __init__(self, specs:dict[str,Any]={}, indexes: [bool] = [], a_string = "GTTACTTGGA", b_string = "GGTTTTTTGC", **kwargs : Any) -> None:
        self.structure = self.init_structure()
        self.specs = specs
        self.indexes = [0]*1728
        # strands from indexes
        self.domains = self.init_domains()
        self.strands = self.indexes2strands(indexes)
        if "max_size" not in self.specs:
            self.specs["max_size"] = 3
        self.complexes = ComplexSet(strands=[s for s,_ in self.strands], complexes=SetSpec(**specs))
        self.evaluation = None
        super().__init__(iterable=self.strands,**kwargs)
        #self.tube = Tube(dict(self.strands), complexes=SetSpec(**specs), name =f"tube_{self.name}")

    def init_domains(self):
        domain_name_lst = ['a', 'b', 'c', 'd', 'e', 'f']
        domain_seq_lst = ["GTTCCAGCACCTTCACT", "TCAGCCGGTGGACTGAG", "GGCAACGTCCTGTTACT", "GAACCCGGGAAAGAC", "TTGAAGAGTAGAGCACA", "CCTAGAGAGGCGCACAT"]
        domains = {}
        for i, domain_name in enumerate(domain_name_lst):
            domains[domain_name] = Domain(domain_seq_lst[i], name=domain_name)
            domains[domain_name + "*"] = ~Domain(domain_seq_lst[i], name=domain_name)
        return domains
    
    def init_structure(self):
        # input/structure_seq/input_seq_L3.csv
        f = open("../../input/structure_seq/input_seq_L3.csv", "r")
        lst = csv.reader(f)
        structure = []
        for l in lst:
            for elem in l:
                structure.append(elem)
        return structure

    def reinit(self):
        self.complexes = ComplexSet(strands=[s for s,_ in self.strands], complexes=SetSpec(**self.specs))
        self.evaluation = None
        self.tube = Tube(dict(self.strands), complexes=SetSpec(**self.specs), name =f"tube")
    
    def indexes2strands(self, indexes, concentration=1e-10):
        tmp_index = 0
        strands = []
        for num, index in enumerate(indexes):
            if index == 1:
                tmp_domain = Domain("", name="s"+str(tmp_index))
                tmp_strand = Strand(str(tmp_domain), name="Strand tmp")
                for domain in self.structure[num].split(" "):
                    tmp_strand = Strand(self.domainlist_to_string(tmp_domain+self.domains[domain]), name="x")
                    # if domain == "a":
                    #     tmp_strand = Strand(self.domainlist_to_string(tmp_domain+self.a), name="x")
                    # elif domain == "b":
                    #     tmp_strand = Strand(self.domainlist_to_string(tmp_domain+self.b), name="x")
                    # elif domain == "a*":
                    #     tmp_strand = Strand(self.domainlist_to_string(tmp_domain+(~(self.a))), name="x")
                    # elif domain == "b*":
                    #     tmp_strand = Strand(self.domainlist_to_string(tmp_domain+(~(self.b))), name="x")
                    tmp_domain = Domain(str(tmp_strand), name="s"+str(tmp_index))                
                self.strands.append((tmp_strand, concentration))
                tmp_index += 1
        return strands
    
    def domainlist_to_string(self, a):
        if isinstance(a,Domain):
            return a.to_string()
        if isinstance(a,DomainList):
            return ''.join([b.to_string() for b in a])

    def __eq__(self, other): 
        if not isinstance(other, L3Individual):
            # don't attempt to compare against unrelated types
            return False
        setA = set(self.strands)
        setB = set(other.strands)
        return setA^setB == set([])

    def eliminateStrand(self):
        nb_strands = [i for i, a in enumerate(self.indexes) if a == 1]
        self.indexes[random.choice(nb_strands)] = 0
        self.reinit()

    def addStrand(self):
        nb_strands = [i for i, a in enumerate(self.indexes) if a == 0]
        self.indexes[random.choice(nb_strands)] = 1
        self.reinit()

@registry.register
class L3GenIndividuals(GenIndividuals):
    def __next__(self):
        return L3Individual()

@registry.register
class L3Evo(Evolution):
    strands_number_domain: DomainLike
    sel_pb: float
    init_pb: float
    mut_pb: float

    def __init__(self, container: Container, budget: int,
            strands_number_domain: DomainLike = (3.,6.), concentration_domain: DomainLike = (1e-10,1e-6),
            sel_pb: float = 0.5, init_pb: float = 0.5, mut_pb: float = 0.2, dimension: int = 2, nupack_specs = {"max_size": 3}, **kwargs):
        self.strands_number_domain = strands_number_domain
        self.concentration_domain = concentration_domain
        self.sel_pb = sel_pb
        self.init_pb = init_pb
        self.mut_pb = mut_pb
        self.nupack_specs = nupack_specs
    
        select_or_initialise = partial(
                tools.sel_or_init,
                sel_fn = tools.sel_random,
                sel_pb = 1,
                init_fn = self.my_init_fn,
                init_pb = 0)
        def vary_fn(ind):
            remove_pb = 0.5 if np.sum(ind.indexes) > strands_number_domain[0] else 0.0
            add_pb = 0.5 if np.sum(ind.indexes) < strands_number_domain[1] else 0.0
            rand_pb = 1 - remove_pb - add_pb
            mut_type = random.choice(3,p=[remove_pb, add_pb, rand_pb])
            if mut_type == 0:
                ind.eliminateStrand()
            elif mut_type == 1:
                ind.addStrand()
            return ind
        
        
        super().__init__(container, budget, dimension=dimension, # type: ignore
                select_or_initialise=select_or_initialise, vary=vary_fn, base_ind_gen=L3GenIndividuals(), **kwargs) # type: ignore

        
    def my_init_fn(self, base_ind):
        # print("called init fn")
        # generate random indexes
        indexes = [0]*1728
        nb_strands = random.choice(range(int(self.strands_number_domain[0]), int(self.strands_number_domain[1])))
        while nb_strands > 0:
            index = random.choice(range(0, 15))
            if indexes[index] == 0:
                indexes[index] = 1
                nb_strands -= 1
        # generate strands from indexes
        base_ind.indexes = indexes
        base_ind.indexes2strands(indexes)
        base_ind.specs = self.nupack_specs
        base_ind.reinit()
        return base_ind

    def _internal_ask(self, base_ind):
        return super()._internal_ask(base_ind)    

def getModel(path="../../saved_model/L3_ave_230530"):
    model = load_model(path)
    return model


# def eval_fn(ind, nupackmodel = Model(material='dna'), fitness_scale = (1e-4)/2,max_strands = 15, complexes_scale=270):
#     """An example evaluation function. It takes an individual as input, and returns the pair ``(fitness, features)``, where ``fitness`` and ``features`` are sequences of scores."""

#     tube_results = nupack.tube_analysis(tubes=[ind.tube], model=nupackmodel)
#     score = 0
#     complexes = tube_results.complexes
#     for t, v in tube_results.tubes.items():
#         #print("tube:",t,v)
#         for c in t.complexes:
#             score += complexes[c].free_energy*v[c]
#     features = ((len(ind.strands)-1)/(max_strands-1), np.tanh((len(ind.tube.complexes)-1)/(complexes_scale-1)))
#     ind.fitness.values = (- np.tanh(score/fitness_scale),)
#     ind.features.values = features
#     #print("DEBUG",ind.fitness, ind.features)
#     return (- np.tanh(score/fitness_scale),), features

def nupack_val(ind, nupackmodel = nupackModel(material='dna'), fitness_scale = (1e-4)/2, complexes_scale=200):
    tube_results = nupack.tube_analysis(tubes=[ind.tube], model=nupackmodel)
    score = 0
    complexes = tube_results.complexes
    for t, v in tube_results.tubes.items():
        #print("tube:",t,v)
        for c in t.complexes:
            score += complexes[c].free_energy*v[c]
    return - np.tanh(score/fitness_scale), np.tanh((len(ind.tube.complexes)-1)/(complexes_scale-1))

def set_eval(ind, averageModel, deviationModel, scale=10.0 ):
    strands = ind.indexes
    # print("nb strands    : ", np.sum(strands))
    # print("index strands : ", [i for i, strand in enumerate(strands) if strand == 1])
    score = 2*math.atan(averageModel.predict([strands], verbose = 0)[0][0]/scale)/math.pi
    fit0 = 2*math.atan(deviationModel.predict([strands], verbose = 0)[0][0]/scale)/math.pi
    energy, nb_comp = nupack_val(ind)
    # fit1 = np.sum(strands)
    fit1 = nb_comp
    # fit1 = energy
    features = (fit0, fit1)
    return (score,), features

def run_qdpy(dirpath="test"):
    # Create container and algorithm. Here we use MAP-Elites, by illuminating a Grid container by evolution.
        
    grid = containers.Grid(
        shape=(16,16), 
        max_items_per_bin=1, 
        fitness_domain=((0.0, 1.),),  
        # features_domain=((0., 1.), (1, 1728))) #軸 deviation, number of strands 
        features_domain=((0., 1.), (0., 1.))) #軸 deviation, free energy
    
    algo = L3Evo(
        grid, 
        budget=10000, 
        batch_size=1000,
        optimisation_task="maximization")
    
    # Create a logger to pretty-print everything and generate output data files
    logger = algorithms.AlgorithmLogger(algo)

    averageModel = getModel('../../saved_model/L3_ave_230613')
    deviationModel = getModel('../../saved_model/L3_dev_230613')
    eval_fn = functools.partial(set_eval,averageModel=averageModel,deviationModel=deviationModel)
    best = algo.optimise(eval_fn)
    print(algo.summary())

    # Plot the results
    logger.final_filename = dirpath + "/qdpy_log_L3_230731.p"
    print(logger.final_filename)
    plots.default_plots_grid(logger)
    print("All results are available in the '%s' pickle file." % logger.final_filename)

def main():
    run_qdpy()
    
if __name__ == "__main__":
    main()
