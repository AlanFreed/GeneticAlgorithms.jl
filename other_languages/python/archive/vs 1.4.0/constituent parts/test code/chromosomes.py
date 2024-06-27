# -*- coding: utf-8 -*-
#!/usr/bin/env python3
import math
import random
import numpy as np
from genes import Gene
"""
    A chromosome is a genetic container of genes.  In this implementation of 
    a genetic algorithm, each chromosome represents a parameter (an unknown 
    value in some model to be parameterized).  Genetic processes (mutation 
    and crossover) adjust the chromosomes and therefore their associated
    parameteric values.  Chromosomes are where genetics and optimization meet.
    A chromosome has an interface of:
    
    Constructor
    
        c = Chromosome(minValue, maxValue, nbrSigFigs)
                minValue    most negative value a parameter can have
                maxValue    most positive value a parameter can have
                nbrSigFigs  seek parameter to this many significant figures
                
    Methods
    
        g = c.pop(l)           returns the gene 'g' held at location 'l'
        c.push(g, l)           assigns gene 'g' to chromosome 'c' at 'l'
        d = c.copy()           returns a deep copy 'd' of chromosome 'c'
        b = c1.isEqualTo(c2)   returns a boolean 'b' regarding truth c1 == c2
        s = c.toString()       returns string 's' describing chromosome 'c'
        c.mutate(probability)  random flip gene expressions at 'probability'
        p = c.decode()         returns the phenotype (parameter) held by 'c'
        c.encode(p)            assigns a phenotype 'p' to chromosome 'c'
        c.genes()              number of genes that comprise this chromosome 
        gA.crossover(gB,pM,pX) Xover between genomes 'gA' and 'gB' with 
                               probabilities of mutation 'pM' and Xover 'pX'
"""


class Chromosome(object):
    
    
    # the two kinds of admissible genes
    recessive = Gene(0)
    dominant = Gene(1)
    
    
    def __init__(self, minValue, maxValue, nbrSigFigs):
        
        if minValue >= maxValue:
            raise RuntimeError(
                "Cannot create a chromosome unless maxValue > minValue.")
        if minValue > 0.0:
            logDecades = math.log10(maxValue / minValue)
        elif maxValue < 0.0:
            logDecades = math.log10(minValue / maxValue)
        elif minValue == 0.0:
            logDecades = math.log10(maxValue)
        else: # maxParameter = 0.0
            logDecades = math.log10(-minValue)
            
        if logDecades > 0.0:
            decades = math.ceil(logDecades)
        else:
            decades = abs(math.floor(logDecades))
        if decades < 1:
            decades = 1
        if decades > 9:
            decades = 9
            
        self.sigFig = nbrSigFigs
            
        if nbrSigFigs == 2:
            bitsPerDecade = 7
        elif nbrSigFigs == 2:
            bitsPerDecade = 10
        elif nbrSigFigs == 3:
            bitsPerDecade = 14
        elif nbrSigFigs == 4:
            bitsPerDecade = 17
        elif nbrSigFigs == 5:
            bitsPerDecade = 20
        elif nbrSigFigs == 6:
            bitsPerDecade = 24
        else:   
            bitsPerDecade = 27
        
        self.bits = bitsPerDecade * decades
            
        self.binary = np.zeros(self.bits, dtype=bool)
        self.gray = np.zeros(self.bits, dtype=bool)
        
        self.chromosome = np.empty(self.bits, dtype=object)
        for i in range(self.bits):
            self.chromosome[i] = Gene()
            
        self.minP = minValue
        self.maxP = maxValue
        self.rangeOfGenes = int(2**(self.bits))
    
        return # new object: this is the constructor
        

    # The encode/decode maps between a real and its Haploid representation.

    #  The algorithm for converting between binary and gray codes assumes the
    #  most significant bit (MSB) is at the left of the code, and associates
    #  with the [1] position in the binary and gray encodings.  The least
    #  significant bit associates with position [High], in other words, e.g.,
    #     code = [1|0|1|1|0|0|1|0] has a MSB of 1 and a LSB of 0. 

    def _binaryToGray(self):
        self.gray[0] = self.binary[0];
        for i in range(1, self.bits):
            self.gray[i] = self.binary[i-1] ^ self.binary[i]
        return # nothing
        
        
    def _grayToBinary(self):
        self.binary[0] = self.gray[0];
        for i in range(1, self.bits):
            self.binary[i] = self.binary[i-1] ^ self.gray[i]
        return # nothing
    
    
    def _binaryToInteger(self): 
        iVal = 0
        iPow = 1
        for i in range(self.bits-1, -1, -1):
            if self.binary[i]:
                iVal = iVal + iPow
            iPow = 2 * iPow
        return iVal
        
    
    def _integerToBinary(self, i):
        atI = i
        bit = self.bits - 1
        while atI > 0:
            if int(atI % 2) == 0:
                self.binary[bit] = False
            else:
                self.binary[bit] = True
            atI = int(atI // 2)
            bit = bit - 1
        # remaining higher-order binary bits are zeros
        for k in range(bit, -1, -1):
            self.binary[k] = False
        return # nothing
        
        
    def _integerToPhenotype(self, i):
        phenotype = (self.minP + (float(i) / float(self.rangeOfGenes))
                     * (self.maxP - self.minP))
        if phenotype < self.minP:
            phenotype = self.minP
        if phenotype > self.maxP:
            phenotype = self.maxP
        return phenotype
        
    
    def _phenotypeToInteger(self, p):
        fraction = (p - self.minP) / (self.maxP - self.minP)
        intVal = int(round(fraction * float(self.rangeOfGenes)))
        if intVal < 0:
            intVal = 0
        maxGene = self.rangeOfGenes
        if intVal > maxGene:
            intVal = maxGene
        return intVal
        
        
    def pop(self, atLocation):
        if atLocation in range(self.bits):
            return self.chromosome[atLocation]
        else:
            raise RuntimeError(
                "The requested gene site to pop from does not exist.")
        
        
    def push(self, gene, toLocation):
        if toLocation in range(self.bits):
            self.chromosome[toLocation].push(gene.pop())
        else:
            raise RuntimeError(
                "The requested gene site to push to does not exist.")
            
            
    def copy(self):
        c = Chromosome(self.minP, self.maxP, self.sigFig)
        for i in range(self.bits):
            c.chromosome[i] = self.chromosome[i].copy()
        return c
        
        
    def isEqualTo(self, c):
        truth = True
        if self.bits == c.bits:
            for i in range(self.bits):
                if not self.chromosome[i].isEqualTo(c.chromosome[i]):
                    truth = False
        else:
            truth = False
        return truth
        
        
    def toString(self):
        l = []
        for i in range(self.bits):
            l.append(self.chromosome[i].toString())
        s = ''.join(l)
        return s
        
        
    def mutate(self, probOfMutation):
        for i in range(self.bits):
            self.chromosome[i].mutate(probOfMutation)
        return # nothing
            
        
    def decode(self):
        for i in range(self.bits):
            if self.chromosome[i].isEqualTo(self.dominant):
                self.gray[i] = True
            else: # recessive
                self.gray[i] = False
        self._grayToBinary()
        j = self._binaryToInteger()
        x = self._integerToPhenotype(j)
        return x
  
  
    def encode(self, phenotype):
        j = self._phenotypeToInteger(phenotype) 
        self._integerToBinary(j)
        self._binaryToGray()
        for i in range(self.bits):
            if self.gray[i] == True:
                self.chromosome[i].push(self.dominant.pop())
            else:
                self.chromosome[i].push(self.recessive.pop())
        
        
    def genes(self):
        n = self.bits
        return n
        
        
"""
    Two chromosomes can reproduce creating an offspring, which is what the
    procedure Crossover mimics.  It accepts two parent chromosomes along with
    probabilities of occurrance for mutation and crossover (a spliting of  
    the chromosome between a random pair of neighboring genes).  The outcome 
    is a child of the parents, and it is returned.  There is a chance that  
    the child will be a clone of one of its parents - a chance that dictated 
    by the assigned probabilities.
"""


def Crossover(parentA, parentB, probOfMutation, probOfCrossover):
    if isinstance(parentA, Chromosome) and isinstance(parentB, Chromosome):
        if parentA.genes() == parentB.genes():
            if random.random() < 0.5:
                child = parentA.copy()
                # left side of chromosome splice belongs to parent A
                if probOfCrossover > random.random():
                    xover = random.randint(3, child.genes()-2)
                    for i in range(xover, child.genes()):
                        child.push(parentB.pop(i), i)
            else:
                child = parentB.copy()
                # left side of chromosome splice belongs to parent B
                if probOfCrossover > random.random():
                    xover = random.randint(3, child.genes()-2)
                    for i in range(xover, child.genes()):
                        child.push(parentA.pop(i), i)
            child.mutate(probOfMutation)
        else:
            raise RuntimeError(
                "Crossover parents must have the same number of genes.")
    else:
        raise RuntimeError("The parents for crossover must be chromosomes.")
    return child