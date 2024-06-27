#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 15:26:33 2017

@author: Prof. Alan D. Freed
"""

import numpy as np
from creatures import Creature


def test():
    vary = np.array([True, False, True])
    fixX = np.array([3.141])
    minX = np.array([-1.0, 1.0])
    optX = np.array([0.1,  5.0])
    maxX = np.array([1.0, 10.0])
    prbM = 0.01
    prbX = 0.85
    sigF = 4
    cre1 = Creature(vary, fixX, minX, maxX, sigF, prbM, prbX)
    cre1.procreate()
    cre2 = Creature(vary, fixX, minX, maxX, sigF, prbM, prbX)
    cre2.alien(optX)
    cre3 = Creature(vary, fixX, minX, maxX, sigF, prbM, prbX)
    cre3.conceive(cre1, cre2)
    print("")
    print("These creatures have {} genes and {} chromosomes."
          .format(cre1.genes(), cre1.chromosomes()))
    print("")
    print("First creature has a gene expression of:")
    print(cre1.toString())
    print("Second creature has a gene expression of:")
    print(cre2.toString())
    print("Their offspring has a gene expression of:")
    print(cre3.toString())
    if cre1.isEqualTo(cre3) or cre2.isEqualTo(cre3):
        print("The child is a clone.")
    else:
        print("The child is unique.")
    print("")
    print("These creatures have parameters:")
    x1 = cre1.getParameters()
    x2 = cre2.getParameters()
    x3 = cre3.getParameters()
    print("For the first:  {}".format(x1))
    print("For the second: {}".format(x2))
    print("For the third:  {}".format(x3))
    print("")
    quality = 0.5
    cre1.setFitness(quality)
    qual = cre1.getFitness()
    print("An assigned fitness of {} was returned as {}."
          .format(quality, qual))


test()
