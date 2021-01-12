#!/usr/bin/env python3
''' Description: Contains viral and population data structures for simulating viral evolution '''

import random
import time
import cProfile
import re

class Virus:
    def __init__(self, genomeSize):
        self.genomeSize = genomeSize
        self.beneficial_d = dict()
        self.neutral_d = dict()
    
    def neutralMutation(self, track_neut):
        base = random.choice([1, 2, 3, 4])
        index = random.randrange(0, self.genomeSize)
        track_neut.append(index)
        self.neutral_d[index] = base
    
    def beneficialMutation(self, track_ben):
        index = random.randrange(0, self.genomeSize)
        track_ben.append(index)
        arr = [1,2,3,4]
        if index in self.neutral_d:
            arr.remove(self.neutral_d[index])
        base = random.choice(arr)
        self.beneficial_d[index] = base

    def __deepcopy__(self, memodict={}):
        copy_object = Virus(self.genomeSize)
        for (key,value) in self.beneficial_d.items():
            copy_object.beneficial_d[key] = value
        for (key,value) in self.neutral_d.items():
            copy_object.neutral_d[key] = value
        return copy_object
    
class Population:
    def __init__(self, r, w, x, probBen, mutRate, genomeSize, numberBeneficial, numberNeutral, initSize, gens, mutations, reference, positiveLoci, neutralLoci):
        self.r = r
        self.w = w
        self.x = x
        self.probBen = probBen
        self.mutRate = mutRate
        self.genomeSize = genomeSize
        self.initSize = initSize
        self.gens = gens
        self.numberBeneficial = numberBeneficial
        self.numberNeutral = numberNeutral
        self.mutations = mutations
        self.reference = reference
        self.positiveLoci = positiveLoci
        self.neutralLoci = neutralLoci