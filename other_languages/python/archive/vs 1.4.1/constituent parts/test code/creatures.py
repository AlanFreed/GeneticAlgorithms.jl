# -*- coding: utf-8 -*-
import numpy as np
from chromosomes import Crossover
from genomes import Genome



"""
    A creature has the genetic material of a genome plus a means for its 
    origin.  After a creature (its object) has been created, it needs to 
    have its genetic material assigned to it.  This can occur in one of 
    three ways: procreation (i.e., random), alien (i.e., divinely assigned) 
    or conceived (i.e., coming from its parents).  The first generation of 
    a colony is procreated with the possible exception that their be an 
    alien (an Adam) assignment.  After that, all creatures are created 
    through conception.  A creature has an interface of:
    
    Constructor
    
        c = Creature(varyPar, fixedPar, minFitPar, maxFitPar, 
                     nbrSigFigs, probMutation, probCrossover)
                varyPar        array of T/F: T if parameter varies, else F
                fixedPar       array of parameters FOR THOSE THAT DO NOT VARY
                               if all parameters vary, then None is to be sent
                minFitPar      array of minimum parameters FOR THOSE THAT VARY
                maxFitPar      array of maximum parameters FOR THOSE THAT VARY
                nbrSigFigs     seek parameters to this many significant figures 
                probMutation   probability of a mutation  event occurring
                probCrossover  probability of a crossover event occurring
                
    Creation is to be immediately followed by a genome assignment scheme
    
        c.procreate()          all genes are assigned by a flipping of a coin
        c.alien(alienPar)      assign best guess at what the parameters may be
        c.conceive(cA, cB)     assign genetic data via crossover and mutation
        
    Methods
    
        n = c.copy()           returns a deep copy 'n' of creature 'c'
        b = c1.isEqualTo(c2)   returns a boolean 'b' regarding truth c1 == c2
        s = c.toString()       returns string 's' describing creature 'c'
        c.genes()              number of genes that comprise creature 'c'
        c.chromosomes()        number of chromosomes in creature 'c'
        c.getParameters()      returns array of all (fixed & varied) parameters
        c.getFitness()         returns the fitness of creature 'c'
        c.setFitness(q)        assign a fitness, or quality 'q' to creature 'c'
"""


class Creature(object):
    
    genome = None
    
    def __init__(self, varyPar, fixedPar, minVaryPar, maxVaryPar, 
                 nbrSigFigs, probOfMutation, probOfCrossover):
                
        self.pars = varyPar.size
        self.varyPs = varyPar
        self.fixedPs = fixedPar
        self.minPs = minVaryPar
        self.maxPs = maxVaryPar
        self.sigFig = nbrSigFigs
        self.probM = probOfMutation
        self.probX = probOfCrossover
            
        return # new object: this is the constructor
        
    
    def _assignParameters(self, assigned):
        self.parameters = np.zeros(self.pars, dtype=float)
        if (self.fixedPs is None) or (self.varyPs is None):
            for i in range(self.pars):
                self.parameters[i] = assigned[i]
        else:
            fixed = -1
            varied = -1
            for i in range(self.pars):
                if self.varyPs[i]:
                    varied = varied + 1
                    self.parameters[i] = assigned[varied]
                else:
                    fixed = fixed + 1
                    self.parameters[i] = self.fixedPs[fixed]
        return # nothing
        
        
    # these are the three constructors for this class
        
        
    def procreate(self):
        self.genome = Genome(self.minPs, self.maxPs, self.sigFig)
        assigned = self.genome.decode()
        self._assignParameters(assigned)
        return # nothing
        
        
    def alien(self, alienPar):
        self.genome = Genome(self.minPs, self.maxPs, self.sigFig)
        self.genome.encode(alienPar)
        self._assignParameters(alienPar)
        return # nothing
        
        
    def conceive(self, parentA, parentB):
        if isinstance(parentA, Creature) and isinstance(parentB, Creature):
            self.genome = Genome(self.minPs, self.maxPs, self.sigFig)
            for i in range(parentA.genome.chromosomes()):
                chromA = parentA.genome.pop(i)
                chromB = parentB.genome.pop(i)
                child = Crossover(chromA, chromB, self.probM, self.probX)
                self.genome.push(child, i)
            assigned = self.genome.decode()
            self._assignParameters(assigned)
        else:
            raise RuntimeError("'parentA' & 'parentB' must be creatures.")
        return # nothing
        
    # the following are its methods
        
    def copy(self):
        c = Creature(self.varyPs, self.fixedPs, self.minPs, self.maxPs, 
                     self.sigFig, self.probM, self.probX)
        c.genome = self.genome.copy()
        assigned = self.genome.decode()
        c._assignParameters(assigned)
        c.quality = self.quality
        return c
        
        
    def isEqualTo(self, c):
        return self.genome.isEqualTo(c.genome)
        
    
    def toString(self):
        return self.genome.toString()
        
        
    def genes(self):
        return self.genome.genes()
        
        
    def chromosomes(self):
        return self.genome.chromosomes()
        
        
    def getParameters(self):
        return self.parameters
        
        
    def getFitness(self):
        return self.quality
        
        
    def setFitness(self, fitness):
        self.quality = fitness
        


"""

class Creature(object):
    
    genome = None
    
    def __init__(self, varyPar, fixedPar, minVaryPar, maxVaryPar, 
                 nbrSigFigs, probOfMutation, probOfCrossover):
                
        self.pars = varyPar.size
        self.varyPs = varyPar
        self.fixedPs = fixedPar
        self.minPs = minVaryPar
        self.maxPs = maxVaryPar
        self.sigFig = nbrSigFigs
        self.probM = probOfMutation
        self.probX = probOfCrossover
            
        return # new object: this is the constructor
        
    
    def _assignParameters(self, assigned):
        self.parameters = np.zeros(self.pars, dtype=float)
        if (self.fixedPs is None) or (self.varyPs is None):
            for i in range(self.pars):
                self.parameters[i] = assigned[i]
        else:
            fixed = -1
            varied = -1
            for i in range(self.pars):
                if self.varyPs[i]:
                    varied = varied + 1
                    self.parameters[i] = assigned[varied]
                else:
                    fixed = fixed + 1
                    self.parameters[i] = self.fixedPs[fixed]
        return # nothing
        
        
    # these are the three constructors for this class
        
        
    def procreate(self):
        self.genome = Genome(self.minPs, self.maxPs, self.sigFig)
        assigned = self.genome.decode()
        self._assignParameters(assigned)
        return # nothing
        
        
    def alien(self, alienPar):
        self.genome = Genome(self.minPs, self.maxPs, self.sigFig)
        self.genome.encode(alienPar)
        self._assignParameters(alienPar)
        return # nothing
        
        
    def conceive(self, parentA, parentB):
        if isinstance(parentA, Creature) and isinstance(parentB, Creature):
            self.genome = Genome(self.minPs, self.maxPs, self.sigFig)
            for i in range(parentA.genome.chromosomes()):
                chromA = parentA.genome.pop(i)
                chromB = parentB.genome.pop(i)
                child = Crossover(chromA, chromB, self.probM, self.probX)
                self.genome.push(child, i)
            assigned = self.genome.decode()
            self._assignParameters(assigned)
        else:
            raise RuntimeError("'parentA' & 'parentB' must be creatures.")
        return # nothing
        
    # the following are its methods
        
    def copy(self):
        c = Creature(self.varyPs, self.fixedPs, self.minPs, self.maxPs, 
                     self.sigFig, self.probM, self.probX)
        c.genome = self.genome.copy()
        assigned = self.genome.decode()
        c._assignParameters(assigned)
        c.quality = self.quality
        return c
        
        
    def isEqualTo(self, c):
        return self.genome.isEqualTo(c.genome)
        
    
    def toString(self):
        return self.genome.toString()
        
        
    def genes(self):
        return self.genome.genes()
        
        
    def chromosomes(self):
        return self.genome.chromosomes()
        
        
    def getParameters(self):
        return self.parameters
        
        
    def getFitness(self):
        return self.quality
        
        
    def setFitness(self, fitness):
        self.quality = fitness
"""