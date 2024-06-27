#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math
import random
import numpy as np
import scipy.stats as st
from species import Specie

"""
Module genAlg.py provides a genetic algorithm for optimization.

Copyright (c) 2018 Alan D. Freed

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""

# Module metadata
__version__ = "1.4.2"
__date__ = "2017-09-20"
__author__ = "Alan D. Freed"
__author_email__ = "adfreed@tamu.edu"


"""
    From the biologic interpretation, the genetic algorithm implemented here
    is a colony of creatures that advances their quality of life (fitness)
    from one generation to the next.  Class GenAlg, defined at the end of
    this file, is a colony or collection of creatures that are sustained
    from one generation to the next.  Mating between creatures occurs through
    a process known as tournament play where the most fit contestant from a
    random selection of contestants is chosen.  Typically, each successive
    generation is more fit than its predicessor, i.e., quality improves.

    It is through class GenAlg that a user applies this genetic algorithm,
    which is why this comment section has been repeated here at the top of
    the file.  This colony of creatures has the following interface:

    Constructor

        c = GenAlg(species, nbrSigFigs, probOfMutations, probOfCrossovers,
                   probOfImmigrants, namePar, varyPar, fixedPar, alienPar,
                   minVaryPar, maxVaryPar)
        where
            species             An instance of class Specie.
                                Here is where the user's model is introduced
                                and the experimental data are assigned.
            nbrSigFigs          The number of significant figures of accuracy
                                sought in the final solution for parameters.
                                There is a strong correlation between this
                                value and how long it takes to get a solution.
                                Values of 4 or 5 are common.
            probOfMutations     The probability of a gene mutation occurring at
                                conception, ie, a gene swapping its expression.
                                Typically this is a low value, e.g., < 0.01.
            probOfCrossovers    The probability of a crossover event occurring
                                at conception, i.e., a chromosome splitting.
                                Typically this is large, e.g., > 0.85.
            probOfImmigrants    The probability of introducing an immigrant
                                into the gene pool after the first generation.
                                (The first generation is all 'immigrants'.)
                                Use as you see fit.  I usually assign it to
                                introduce one immigrant per generation or so.
            namePar             A string associated with each model parameter.
                                This makes the report much more readable.
            varyPar             Specifies if a particular parameter is to be
                                varied (True) by the optimizer or held fixed
                                (False).  If all parameters vary you can also
                                send None.
            fixedPar            An array of the fixed parameters.  If all
                                parameters vary then send None.
            alienPar            If the user has a 'best guess', then this is
                                where those parameters are input.  If the user
                                does not, then send None.
            minVaryPar          For those parameters that are to be varied by
                                the genetic algorithm, these values represent
                                their respective lower bounds.
            maxVaryPar          For those parameters that are to be varied by
                                the genetic algorithm, these values represent
                                their respective upper bounds.

        The above arguments can accept an input value of None:
            varyPar, fixedPar, alienPar

        IMPORTANT!!!!!
        For the parameter arrays, index 0 is reserved for selecting the best
        objective function.  Specifically, a dimensionless error is computed
            error[i] = (|model[i] - experiment[i]| / sdtDev(experiment))**p
        where i indexes over the set of all experimental data being considered
        and parameter p (> 0) is a characteristic of the objective function.
        It is stored at the 0 index in the various arrays for parameter data.
        Once the error vector is populated the mean error is calculated via
            meanError = mean(error)**(1/p)
        thereby allowing for a quality value to be assigned as
            quality = 1 / meanError
        where the creature whose quality is the greatest is deemed most fit.
        The exponent 'p' is treated as another parameter to be optimized.  It
        can be fixed or allowed to vary.  Setting fixedPar[0] = 2 would impose
        a root mean squared error (RMSE).  Setting fixedPar[0] = 1 imposes a
        linear error.  Or setting varyPar[0] = True with minVaryPar[0] = 0.01
        and maxVaryPar[0] = 100 would seek out an optimal objective function
        in the range of p within [0.01, 100] for the specific set of data that
        are being considered.  I.e., it provides an optimum objective function.

    Methods

        p = c.bestCreatureParam()
            p       The full set of model parameters associated with the most
                    fit or elite creature in the colony at this generation.

        q = c.bestCreatureFitness()
            q       The fitness parameter associated with the most fit or
                    elite creature in the colony.  It is a measure of quality,
                    specifically:
                        quality[i] = 1 / |modelResponse[i] - experiment[i]|
                        quality    = sum over all quality[i]

        n = c.population()
            n       The number of creatures that comprise a single generation.

        g = c.atGeneration()
            g       The number of generations that have been in existance.

        n = c.generationsToConvergence()
            n       The expected number of generations needed for the genetic
                    algorithm to converge upon a solution.

        c.advanceToNextGeneration()
                    Advances the genetic algorithm to its next generation.
                    Giving the user access to this mechanism allows him/her
                    control over how they choose to proceed, and what
                    information they may wish to capture with each successive
                    interation in solution refinement.

        s = c.report()
                    This method returns a string 's' that comprises a report
                    on the status of the optimizer at its current generation.
                    It can/should be called immediately after a colony object
                    is created and right after each advanceToNextGeneration
                    call is made.  A wealth of statistics are reported here,
                    specifically: optimal value, arithmatic mean, median,
                    standard deviation, skewness and excess kurtosis.  This
                    will provide the user with a means to assess viability
                    of the returned solution, i.e., to pass judgment.
"""


#   Proceeding with the interfaces and implementations of the various classes
#   that assemble into a genetic algorithm useful for parameter estimation.


"""
    Genes are the lowest level containers of genetic information.  Here they
    are considered to contain haploid gene expressions: dominant (assigned 1)
    and recessive (assigned 0).  A gene has the following interface:

    Constructors

        g = Gene()     creates a new gene 'g' with a  random   gene expression
        g = Gene(0)    creates a new gene 'g' with a recessive gene expression
        g = Gene(1)    creates a new gene 'g' with a dominant  gene expression

    Methods

        e = g.pop()            returns the gene expression 'e' held in 'g'
        g.push(e)              assigns the gene expression 'e' to gene 'g'
        c = g.copy()           returns a deep copy 'c' of gene 'g'
        b = g1.isEqualTo(g2)   returns a boolean 'b' regarding truth g1 == g2
        s = g.toString()       returns a string representation 's' of gene 'g'
        g.mutate(probability)  random flip of gene expression at 'probability'
"""


class Gene(object):

    # base gene type: 0 is recessive, 1 is dominant
    haploid = (0, 1)

    def __init__(self, geneExpression=None):
        if geneExpression is None:
            self.allele = random.choice(self.haploid)
        elif geneExpression in self.haploid:
            self.allele = geneExpression
        else:
            self.allele = random.choice(self.haploid)
        return  # object: this is the constructor

    def pop(self):
        return self.allele

    def push(self, geneExpression):
        if geneExpression in self.haploid:
            self.allele = geneExpression
        else:
            raise RuntimeError(
                     "A pushed gene expression must be haploid, i.e., 0 or 1.")

    def copy(self):
        g = Gene(self.allele)
        return g

    def isEqualTo(self, gene):
        if gene.allele == self.allele:
            return True
        else:
            return False

    def toString(self):
        if self.allele == 0:
            return '0'
        else:
            return '1'

    # the probability of mutation must lie within (0,1)
    def mutate(self, probabilityOfMutation):
        if probabilityOfMutation > random.random():
            if self.allele == 0:
                self.allele = 1
            else:
                self.allele = 0
        return  # nothing


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
        else:  # maxParameter = 0.0
            logDecades = math.log10(-minValue)

        if logDecades > 0.0:
            decades = int(math.ceil(logDecades))
        else:
            decades = abs(int(math.floor(logDecades)))
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

        return  # new object: this is the constructor

    # The encode/decode maps between a real and its Haploid representation.

    #  The algorithm for converting between binary and gray codes assumes the
    #  most significant bit (MSB) is at the left of the code, and associates
    #  with the [1] position in the binary and gray encodings.  The least
    #  significant bit associates with position [High], in other words, e.g.,
    #     code = [1|0|1|1|0|0|1|0] has a MSB of 1 and a LSB of 0.

    def _binaryToGray(self):
        self.gray[0] = self.binary[0]
        for i in range(1, self.bits):
            self.gray[i] = self.binary[i-1] ^ self.binary[i]
        return  # nothing

    def _grayToBinary(self):
        self.binary[0] = self.gray[0]
        for i in range(1, self.bits):
            self.binary[i] = self.binary[i-1] ^ self.gray[i]
        return  # nothing

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
        return  # nothing

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
        allele = []
        for i in range(self.bits):
            allele.append(self.chromosome[i].toString())
        s = ''.join(allele)
        return s

    def mutate(self, probOfMutation):
        for i in range(self.bits):
            self.chromosome[i].mutate(probOfMutation)
        return  # nothing

    def decode(self):
        for i in range(self.bits):
            if self.chromosome[i].isEqualTo(self.dominant):
                self.gray[i] = True
            else:  # recessive
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
            if self.gray[i]:
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

        return  # new object: this is the constructor

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
                    self.genome[toLocation].push(chromosome.pop(i), i)
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
            s = s + self.genome[i].toString()
            if i < self.vecs-1:
                s = s + '\n'  # line feed
        return s

    def mutate(self, probOfMutation):
        for i in range(self.vecs):
            self.genome[i].mutate(probOfMutation)
        return  # nothing

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

        return  # new object: this is the constructor

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
        return  # nothing

    # these are the three constructors for this class

    def procreate(self):
        self.genome = Genome(self.minPs, self.maxPs, self.sigFig)
        assigned = self.genome.decode()
        self._assignParameters(assigned)
        return  # nothing

    def alien(self, alienPar):
        self.genome = Genome(self.minPs, self.maxPs, self.sigFig)
        self.genome.encode(alienPar)
        self._assignParameters(alienPar)
        return  # nothing

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
        return  # nothing

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
    From the biologic interpretation, the genetic algorithm implemented here
    is a colony of creatures that advances their quality of life (ie, fitness)
    from one generation to the next.  Class GenAlg, defined at the end of
    this file, is a colony or collection of creatures that are sustained
    from one generation to the next.  Mating between creatures occurs through
    a process known as tournament play where the most fit contestant from a
    random selection of contestants is chosen.  Typically, each successive
    generation is more fit than its predicessor, i.e., quality improves.

    It is through class GenAlg that a user applies this genetic algorithm,
    which is why this comment section has been repeated here at the top of
    the file.  This colony of creatures has the following interface:

    Constructor

        c = GenAlg(species, nbrSigFigs, probOfMutations, probOfCrossovers,
                   probOfImmigrants, namePar, varyPar, fixedPar, alienPar,
                   minVaryPar, maxVaryPar)
        where
            species             An instance of class Species (at top of file).
                                Here is where the user's model is introduced
                                and the experimental data are assigned.
            nbrSigFigs          The number of significant figures of accuracy
                                sought in the final solution for parameters.
                                There is a strong correlation between this
                                value and how long it takes to get a solution.
                                Values of 4 or 5 are common.
            probOfMutations     The probability of a gene mutation occurring at
                                conception, ie, a gene swapping its expression.
                                Typically this is a low value, e.g., < 0.01.
            probOfCrossovers    The probability of a crossover event occurring
                                at conception, i.e., a chromosome splitting.
                                Typically this is large, e.g., > 0.85.
            probOfImmigrants    The probability of introducing an immigrant
                                into the gene pool after the first generation.
                                (The first generation is all 'immigrants'.)
                                Use as you see fit.  I usually assign it to
                                introduce one immigrant per generation or so.
            namePar             A string associated with each model parameter.
                                This makes the report much more readable.
            varyPar             Specifies if a particular parameter is to be
                                varied (T) by the optimizer or held fixed (F).
                                If all parameters vary you can also send None.
            fixedPar            An array of the fixed parameters.  If all
                                parameters vary then send None.
            alienPar            If the user has a 'best guess', then this is
                                where those parameters are input.  If the user
                                does not, then send None.
            minVaryPar          For those parameters that are to be varied by
                                the genetic algorithm, these values represent
                                their respective lower bounds.
            maxVaryPar          For those parameters that are to be varied by
                                the genetic algorithm, these values represent
                                their respective upper bounds.

        The above arguments can accept an input value of None:
            varyPar, fixedPar, alienPar

        IMPORTANT!!!!!
        For the parameter arrays, index 0 is reserved for selecting the best
        objective function.  Specifically, a dimensionless error is computed
            error[i] = (|model[i] - experiment[i]| / sdtDev(experiment))**p
        where i indexes over the set of all experimental data being considered
        and parameter p (> 0) is a characteristic of the objective function.
        It is stored at the 0 index in the various arrays for parameter data.
        Once the error vector is populated the mean error is calculated via
            meanError = mean(error)**(1/p)
        thereby allowing for a quality value to be assigned as
            quality = 1 / meanError
        where the creature whose quality is the greatest is deemed most fit.
        The exponent 'p' is treated as another parameter to be optimized.  It
        can be fixed or allowed to vary.  Setting fixedPar[0] = 2 would impose
        a root mean squared error (RMSE).  Setting fixedPar[0] = 1 imposes a
        linear error.  Or setting varyPar[0] = True with minVaryPar[0] = 0.01
        and maxVaryPar[0] = 100 would seek out an optimal objective function
        in the range of p within [0.01, 100] for the specific set of data that
        are being considered.  I.e., it provides an optimum objective function.

    Methods

        p = c.bestCreatureParam()
            p       The full set of model parameters associated with the most
                    fit or elite creature in the colony at this generation.


        q = c.bestCreatureFitness()
            q       The fitness parameter associated with the most fit or
                    elite creature in the colony.  It is a measure of quality,
                    specifically:
                        quality[i] = 1 / |modelResponse[i] - experiment[i]|
                        quality    = sum over all quality[i]


        n = c.population()
            n       The number of creatures that comprise a single generation.


        g = c.atGeneration()
            g       The number of generations that have been in existance.


        n = c.generationsToConvergence()
            n       The expected number of generations needed for the genetic
                    algorithm to converge upon a solution.


        c.advanceToNextGeneration()
                    Advances the genetic algorithm to its next generation.
                    Giving the user access to this mechanism allows him/her
                    control over how they choose to proceed, and what
                    information they may wish to capture with each successive
                    interation in solution refinement.

        s = c.report()
                    This method returns a string 's' that comprises a report
                    on the status of the optimizer at its current generation.
                    It can/should be called immediately after a colony object
                    is created and right after each advanceToNextGeneration
                    call is made.  A wealth of statistics are reported here,
                    specifically: optimal value, arithmatic mean, median,
                    standard deviation, skewness and excess kurtosis.  This
                    will provide the user with a means to assess viability
                    of the returned solution, i.e., to pass judgment.
"""


class GenAlg(object):

    alphabet = 2     # haploid genes have two expressions: dominant & recessive

    def __init__(self, species, nbrSigFigs, probOfMutations, probOfCrossovers,
                 probOfImmigrants, namePar, varyPar, fixedPar, alienPar,
                 minVaryPar, maxVaryPar):

        if isinstance(species, Specie):
            self.species = species
        else:
            raise RuntimeError(
                        "Supplied 'species' does not inherit class 'Species'.")

        # construct population variances for the experimental response data
        self.stdDev = np.empty(self.species.nExp, dtype=object)
        for i in range(self.species.nExp):
            self.stdDev[i] = np.zeros(self.species.nRsp[i], dtype=float)
            for j in range(self.species.nRsp[i]):
                self.stdDev[i][j] = np.std(self.species.expR[i][j, :])

        if nbrSigFigs < 1:
            self.sigFig = 1
        elif nbrSigFigs > 7:
            self.sigFig = 7
        else:
            self.sigFig = nbrSigFigs

        if probOfMutations < 1.0e-6:
            self.probM = 1.0e-6
        elif probOfMutations > 0.999999:
            self.probM = 0.999999
        else:
            self.probM = probOfMutations

        if probOfCrossovers < 1.0e-6:
            self.probX = 1.0e-6
        elif probOfCrossovers > 0.999999:
            self.probX = 0.999999
        else:
            self.probX = probOfCrossovers

        # immigrants are procreated
        if probOfImmigrants < 1.0e-6:
            self.probI = 1.0e-6
        elif probOfImmigrants > 0.999999:
            self.probI = 0.999999
        else:
            self.probI = probOfImmigrants

        if isinstance(namePar, np.ndarray):
            dimP = namePar.size
            self.namePs = namePar
        else:
            raise RuntimeError("'namePar' must be a numpy array.")

        if varyPar is None:
            self.varyPs = None
            fixedPar = None
        else:
            if (isinstance(varyPar, np.ndarray) and
               (namePar.size == varyPar.size)):
                self.varyPs = varyPar
            else:
                raise RuntimeError("'varyPar' must have length of 'namePar'.")

        fixP = 0
        if fixedPar is None:
            self.fixedPs = None
            varP = dimP
        else:
            varP = 0
            for i in range(dimP):
                if self.varyPs[i]:
                    varP = varP + 1
                else:
                    fixP = fixP + 1
            if (isinstance(fixedPar, np.ndarray) and
               (fixedPar.size == fixP)):
                self.fixedPs = fixedPar
            else:
                raise RuntimeError(
                    "'fixedPar' must be a numpy array of length {}."
                    .format(fixP))

        if alienPar is not None:
            if not (isinstance(alienPar, np.ndarray) and
               (alienPar.size == varP)):
                    raise RuntimeError(
                               "'alienPar' must be a numpy array of length {}."
                               .format(varP))

        if (isinstance(minVaryPar, np.ndarray) and
           (minVaryPar.size == varP)):
            self.minPs = minVaryPar
        else:
            raise RuntimeError("'minVaryPar' must be numpy array of length {}."
                               .format(varP))

        if (isinstance(maxVaryPar, np.ndarray) and
           (maxVaryPar.size == varP)):
            self.maxPs = maxVaryPar
        else:
            raise RuntimeError("'maxVaryPar' must be numpy array of length {}."
                               .format(varP))

        # introduce an alien (or an Adam) into the population
        self.elite = Creature(self.varyPs, self.fixedPs, self.minPs,
                              self.maxPs, self.sigFig, self.probM, self.probX)
        if alienPar is None:
            self.elite.procreate()
        else:
            self.elite.alien(alienPar)
        self._evaluate(self.elite)

        # determine the population size
        dimSchemata = self.sigFig
        # algorithm of D. Goldberg (2002) for estimating population size
        self.popSize = (int(math.ceil(self.alphabet**dimSchemata
                        * dimSchemata * math.log(self.alphabet)
                        + math.log(self.elite.genes()))))

        # create matrix that holds all model parameters: used for statistics
        self.Ps = np.zeros(shape=(varP, self.popSize), dtype=float)
        p = self.elite.getParameters()
        v = -1
        for i in range(dimP):
            if (self.fixedPs is None) or self.varyPs[i]:
                v = v + 1
                self.Ps[v, 0] = p[i]

        # create the population
        self.adults = np.empty(self.popSize, dtype=object)
        self.fitness = np.zeros(self.popSize, dtype=float)
        self.adults[0] = self.elite
        self.fitness[0] = self.elite.getFitness()
        for j in range(1, self.popSize):
            self.adults[j] = Creature(self.varyPs, self.fixedPs, self.minPs,
                                      self.maxPs, self.sigFig, self.probM,
                                      self.probX)
            self.adults[j].procreate()
            self._evaluate(self.adults[j])
            self.fitness[j] = self.adults[j].getFitness()
            p = self.adults[j].getParameters()
            v = -1
            for i in range(dimP):
                if (self.fixedPs is None) or (self.varyPs[i] is True):
                    v = v + 1
                    self.Ps[v, j] = p[i]
        self.generation = 1

        # select the most fit in the population: the elite is at location [0]
        quality = self.elite.getFitness()
        atAdult = 0
        for i in range(1, self.popSize):
            if self.adults[i].getFitness() > quality:
                quality = self.adults[i].getFitness()
                atAdult = i
        if atAdult > 0:
            self.adults[0] = self.adults[atAdult]
            self.adults[atAdult] = self.elite
            self.elite = self.adults[0]
            temp = self.fitness[0]
            self.fitness[0] = self.fitness[atAdult]
            self.fitness[atAdult] = temp

        # data required for tournament play
        self.combatants = int(self.popSize // 50)
        if self.combatants < 3:
            self.combatants = 3
        self.contestants = np.empty(self.combatants, dtype=object)
        self.children = np.empty(self.popSize, dtype=object)

        return  # new object: this is the constructor

    #  select parents
    def _tournamentPlay(self):
        if self.probI > random.random():
            mostFit = Creature(self.varyPs, self.fixedPs, self.minPs,
                               self.maxPs, self.sigFig, self.probM, self.probX)
            mostFit.procreate()
            self._evaluate(mostFit)
        else:
            for i in range(self.combatants):
                combatant = random.randint(0, self.popSize-1)
                self.contestants[i] = self.adults[combatant]
            mostFit = self.contestants[0]
            for i in range(1, self.combatants):
                if self.contestants[i].getFitness() > mostFit.getFitness():
                    mostFit = self.contestants[i]
        return mostFit

    # conceive the children for the next generation
    def _mate(self):
        # elite creature from current generation lives into next generation
        self.children[0] = self.elite
        self.fitness[0] = self.elite.getFitness()
        p = self.children[0].getParameters()
        v = -1
        for i in range(p.size):
            if (self.fixedPs is None) or self.varyPs[i]:
                v = v + 1
                self.Ps[v, 0] = p[i]
        for j in range(1, self.popSize):
            parentA = self._tournamentPlay()
            parentB = self._tournamentPlay()
            count = 0
            while parentA.isEqualTo(parentB):
                parentB = self._tournamentPlay()
                count = count + 1
                if count > int(self.popSize // 2):
                    print("The population is filling up with clones.")
                    print("Consider terminating your run.")
            self.children[j] = Creature(self.varyPs, self.fixedPs, self.minPs,
                                        self.maxPs, self.sigFig, self.probM,
                                        self.probX)
            self.children[j].conceive(parentA, parentB)
            self._evaluate(self.children[j])
            self.fitness[j] = self.children[j].getFitness()
            p = self.children[j].getParameters()
            v = -1
            for i in range(p.size):
                if (self.fixedPs is None) or self.varyPs[i]:
                    v = v + 1
                    self.Ps[v, j] = p[i]
        # select most fit in the new population: the elite is at location [0]
        quality = self.fitness[0]
        atChild = 0
        for i in range(1, self.popSize):
            if self.fitness[i] > quality:
                quality = self.fitness[i]
                atChild = i
        if atChild > 0:
            child = self.children[0]
            self.children[0] = self.children[atChild]
            self.children[atChild] = child
            temp = self.fitness[0]
            self.fitness[0] = self.fitness[atChild]
            self.fitness[atChild] = temp
        return  # nothing

    #  determine the fitness of the creature
    def _evaluate(self, creature):
        par = creature.getParameters()
        self.species.model(par)
        if ((par[0] < np.finfo(float).eps) and (par[0] >= 0)):
            par[0] = np.finfo(float).eps
        if ((par[0] > -np.finfo(float).eps) and (par[0] < 0)):
            par[0] = -np.finfo(float).eps
        error = np.empty(self.species.nExp, dtype=object)
        for i in range(self.species.nExp):
            error[i] = np.zeros(shape=(self.species.nRsp[i],
                                       self.species.nPts[i]), dtype=float)
        n = 0
        quality = 0.0
        for i in range(self.species.nExp):
            for j in range(self.species.nRsp[i]):
                n = n + 1
                for k in range(self.species.nPts[i]):
                    error[i][j, k] = ((abs(self.species.modR[i][j, k]
                                           - self.species.expR[i][j, k])
                                      / self.stdDev[i][j])**par[0])
                meanErr = np.mean(error[i][j, :])**(1.0 / par[0])
                quality = quality + 1.0 / meanErr
        quality = quality / n
        creature.setFitness(quality)
        return  # nothing

    def bestCreatureParam(self):
        return self.elite.getParameters()

    def bestCreatureFitness(self):
        return self.elite.getFitness()

    def population(self):
        return self.popSize

    def atGeneration(self):
        return self.generation

    # algorithm of D. Goldberg (2002) for estimating convergence
    def generationsToConvergence(self):
        convergeAt = int(math.ceil(math.sqrt(self.elite.genes())
                         * math.log(self.popSize) / math.log(self.combatants)))
        return convergeAt

    def advanceToNextGeneration(self):
        self.generation = self.generation + 1
        self._mate()
        for i in range(self.popSize):
            self.adults[i] = self.children[i]
        self.elite = self.adults[0]
        return  # nothing

    def report(self):
        s = '\n' + 'Statistics for generation {}'.format(self.generation)
        s = s + ' with a population size of {}.'.format(self.popSize)
        s = s + '\n'
        s = s + 'Optimal value and population statistics for fitness:'
        s = s + '\n'
        fitness = self.elite.getFitness()
        if fitness >= 0.0:
            s = s + '   optimal value    {:.4E}'.format(fitness)
        else:
            s = s + '   optimal value   {:.4E}'.format(fitness)
        s = s + '\n'
        mean = np.mean(self.fitness)
        if mean >= 0.0:
            s = s + '   arithmatic mean  {:.4E}'.format(mean)
        else:
            s = s + '   arithmatic mean {:.4E}'.format(mean)
        s = s + '\n'
        median = np.median(self.fitness)
        if median >= 0.0:
            s = s + '   median           {:.4E}'.format(median)
        else:
            s = s + '   median          {:.4E}'.format(median)
        s = s + '\n'
        std = np.std(self.fitness)
        if std >= 0.0:
            s = s + '   std deviation    {:.4E}'.format(std)
        else:
            s = s + '   std deviation   {:.4E}'.format(std)
        s = s + '\n'
        skew = st.skew(self.fitness)
        if skew >= 0.0:
            s = s + '   skewness         {:.4E}'.format(skew)
        else:
            s = s + '   skewness        {:.4E}'.format(skew)
        s = s + '\n'
        kurtosis = st.kurtosis(self.fitness) - 3
        if kurtosis >= 0.0:
            s = s + ('   excess kurtosis  {:.4E}'.format(kurtosis))
        else:
            s = s + ('   excess kurtosis {:.4E}'.format(kurtosis))
        s = s + '\n'
        p = self.elite.getParameters()
        v = -1
        for i in range(p.size):
            if (self.fixedPs is None) or self.varyPs[i]:
                v = v + 1
                s = s + 'Optimal value and population statistics '
                s = s + 'for parameter '
                s = s + self.namePs[i]
                s = s + ':'
                s = s + '\n'
                param = p[i]
                if param >= 0.0:
                    s = s + '   optimal value    {:.4E}'.format(param)
                else:
                    s = s + '   optimal value   {:.4E}'.format(param)
                s = s + '\n'
                mean = np.mean(self.Ps[v, :])
                if mean >= 0.0:
                    s = s + '   arithmatic mean  {:.4E}'.format(mean)
                else:
                    s = s + '   arithmatic mean {:.4E}'.format(mean)
                s = s + '\n'
                median = np.median(self.Ps[v, :])
                if median >= 0.0:
                    s = s + '   median           {:.4E}'.format(median)
                else:
                    s = s + '   median          {:.4E}'.format(median)
                s = s + '\n'
                std = np.std(self.Ps[v, :])
                if std >= 0.0:
                    s = s + '   std deviation    {:.4E}'.format(std)
                else:
                    s = s + '   std deviation   {:.4E}'.format(std)
                s = s + '\n'
                skew = st.skew(self.Ps[v, :])
                if skew >= 0.0:
                    s = s + '   skewness         {:.4E}'.format(skew)
                else:
                    s = s + '   skewness        {:.4E}'.format(skew)
                s = s + '\n'
                kurtosis = st.kurtosis(self.Ps[v, :]) - 3
                if kurtosis >= 0.0:
                    s = s + ('   excess kurtosis  {:.4E}'.format(kurtosis))
                else:
                    s = s + ('   excess kurtosis {:.4E}'.format(kurtosis))
                s = s + '\n'
        s = s + 'Genetic genome for the most fit creature:'
        s = s + '\n'
        s = s + self.elite.toString()
        s = s + '\n'
        if self.fixedPs is not None:
            s = s + 'Parameters assigned fixed values:'
            s = s + '\n'
            f = -1
            for i in range(p.size):
                if not self.varyPs[i]:
                    f = f + 1
                    s = s + '   '
                    s = s + self.namePs[i]
                    s = s + ':  {:.4E}'.format(self.fixedPs[f])
                    s = s + '\n'
        return s
