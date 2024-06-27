#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import random


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
