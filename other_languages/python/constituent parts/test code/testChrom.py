#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from chromosomes import Chromosome
from chromosomes import Crossover
"""
Created on Mon Jul 10 08:46:11 2017

@author: Prof. Alan D. Freed
"""


def test():
    minX = 1.0
    maxX = 5.0
    prbM = 0.01
    prbX = 0.85
    sigF = 4
    crm1 = Chromosome(minX, maxX, sigF)
    crm2 = Chromosome(minX, maxX, sigF)
    crm3 = Crossover(crm1, crm2, prbM, prbX)
    print("These chromosomes are made of {} genes.".format(crm1.genes()))
    print("First  chromosome has gene expression {}.".format(crm1.toString()))
    print("Second chromosome has gene expression {}.".format(crm2.toString()))
    print("Their offspring has a gene expression {}.".format(crm3.toString()))
    if crm1.isEqualTo(crm3) or crm2.isEqualTo(crm3):
        print("   The child is a clone.")
    else:
        print("   The child is unique.")
    print("These chromosomes decode as:")
    x1 = crm1.decode()
    x2 = crm2.decode()
    x3 = crm3.decode()
    print("   For the first:  {}".format(x1))
    print("   For the second: {}".format(x2))
    print("   For the third:  {}".format(x3))
    print("And then encode as:")
    crm1.encode(x1)
    crm2.encode(x2)
    crm3.encode(x3)
    print("   For the first:  {}".format(crm1.toString()))
    print("   For the second: {}".format(crm2.toString()))
    print("   For the third:  {}".format(crm3.toString()))
    print("The 5th gene of this sequence of chromosomes is: {}, {}, {}."
          .format(crm1.pop(5).pop(), crm2.pop(5).pop(), crm3.pop(5).pop()))
    c = crm3.copy()
    print("A copy of the third chromosome is {}.".format(c.toString()))
    c.push(Chromosome.dominant, 5)
    print("Its fifth gene is made dominant   {}.".format(c.toString()))


test()
