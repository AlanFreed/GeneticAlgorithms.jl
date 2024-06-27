# -*- coding: utf-8 -*-
import numpy as np
from chromosomes import Chromosome
from chromosomes import Crossover
"""
    A genome is a genetic container of chromosomes.  Here is where all the 
    genetic information that makes a creature a creature exists.  A genome 
    has an interface of:
    
    Constructor
    
        c = Genome(minValues, maxValues, nbrSigFigs)
                minValues      array: most negative values parameters can have
                maxValues      array: most positive values parameters can have
                nbrSigFigs     seek parameters to this many significant figures 
                
    Methods
    
        c = g.pop(l)           returns the chromosome 'c' held at location 'l'
        g.push(c, l)           assigns chromosome 'c' to genome 'g' at 'l'
        c = g.copy()           returns a deep copy 'c' of genome 'g'
        b = g1.isEqualTo(g2)   returns a boolean 'b' regarding truth g1 == g2
        s = g.toString()       returns string 's' describing genome 'g'
        c.mutate(probability)  random flip gene expressions at 'probability'
        p = g.decode()         returns the phenotypes (parameters) held by 'g'
        g.encode(arrayP)       assigns phenotypes 'arrayP' to genome 'g'
        g.genes()              number of genes that comprise genome 'g'
        g.chromosomes()        number of chromosomes that are in the genome 'g'
"""
    

class Genome(object):
    
    def __init__(self, minValues, maxValues, nbrSigFigs):
        
        self.sigFig = nbrSigFigs
        self.vecs = maxValues.size
        self.minPs = minValues
        self.maxPs = maxValues
                    
        self.genome = np.empty(self.vecs, dtype=object)
        for i in range(self.vecs):
            minP = self.minPs[i]
            maxP = self.maxPs[i]
            self.genome[i] = Chromosome(minP, maxP, nbrSigFigs)
            
        return # new object: this is the constructor
        
        
    def pop(self, atLocation):
        if atLocation in range(self.vecs):
            return self.genome[atLocation]
        else:
            raise RuntimeError(
                "The requested chromosome site to pop from does not exist.")
        
        
    def push(self, chromosome, toLocation):
        if toLocation in range(self.vecs):
            if self.genome[toLocation].genes() == chromosome.genes():
                for i in range(chromosome.genes()):
                    self.genome[toLocation].push(chromosome.pop(i),i)
            else:
                raise RuntimeError(
                    "Chromosome has wrong length for pushing into genome.")
        else:
            raise RuntimeError(
                "The requested chromosome site to push to does not exist.")
        
        
    def copy(self):
        g = Genome(self.minPs, self.maxPs, self.sigFig)
        for i in range(self.vecs):
            g.genome[i] = self.genome[i].copy()
        return g
        
        
    def isEqualTo(self, g):
        truth = True
        if self.chromosomes() == g.chromosomes():
            for i in range(self.vecs):
                if not self.genome[i].isEqualTo(g.genome[i]):
                    truth = False
        else:
            truth = False
        return truth
                
        
    def toString(self):
        s = ''
        for i in range(self.vecs):
            s = s + ('{}: '.format(i))
            s = s + self.genome[i].toString();
            if i < self.vecs-1:
                s = s + '\n'  # line feed
        return s
        
        
    def mutate(self, probOfMutation):
        for i in range(self.vecs):
            self.genome[i].mutate(probOfMutation)
        return # nothing
            
        
    def decode(self):
        phenotypes = np.zeros(self.vecs, dtype=float)
        for i in range(self.vecs):
            phenotypes[i] = self.genome[i].decode()
        return phenotypes

  
    def encode(self, phenotypes):
        if phenotypes.size == self.vecs:
            for i in range(self.vecs):
                self.genome[i].encode(phenotypes[i])
        else:
            raise RuntimeError("Phenotypes for encoding had the wrong length.")
        
    def genes(self):
        genes = 0
        for i in range(self.vecs):
            genes = genes + self.genome[i].genes()
        return genes
        
        
    def chromosomes(self):
        c = self.vecs
        return c
     
    def crossover(self, parentB, probOfMutation, probOfCrossover):
        if isinstance(parentB, Genome):
            if self.chromosomes() == parentB.chromosomes():
                child = self.copy()
                for i in range(child.chromosomes()):
                    chromosome = Crossover(self.pop(i), parentB.pop(i), 
                                           probOfMutation, probOfCrossover)
                    child.push(chromosome, i)
            else:
                raise RuntimeError(
                "Crossover parents must have the same number of chromosomes.")
        else:
            raise RuntimeError(
                "The parents for crossover must both be genomes.")
        return child