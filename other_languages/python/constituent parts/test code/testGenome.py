#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 13:06:01 2017

@author: Prof. Alan D. Freed
"""

import numpy as np
from genomes import Genome


def test():
    minX = np.array([0.0, -1.0, -5.0])
    maxX = np.array([1.0, 1.0, -1.0])
    prbM = 0.01
    prbX = 0.8
    sigF = 4
    gnm1 = Genome(minX, maxX, sigF)
    gnm2 = Genome(minX, maxX, sigF)
    gnm3 = gnm1.crossover(gnm2, prbM, prbX)
    print("")
    print("These genomes have {} genes and {} chromosomes."
          .format(gnm1.genes(), gnm1.chromosomes()))
    print("")
    print("First genome has a gene expression of:")
    print(gnm1.toString())
    print("Second genome has a gene expression of:")
    print(gnm2.toString())
    print("Their offspring has a gene expression of:")
    print(gnm3.toString())
    if gnm1.isEqualTo(gnm3) or gnm2.isEqualTo(gnm3):
        print("The child is a clone.")
    else:
        print("The child is unique.")
    print("")
    print("These genomes decode as:")
    x1 = gnm1.decode()
    x2 = gnm2.decode()
    x3 = gnm3.decode()
    print("For the first:  {}".format(x1))
    print("For the second: {}".format(x2))
    print("For the third:  {}".format(x3))
    print("And then encode as:")
    print("")
    gnm1.encode(x1)
    gnm2.encode(x2)
    gnm3.encode(x3)
    print("First genome has a gene expression of:")
    print(gnm1.toString())
    print("Second genome has a gene expression of:")
    print(gnm2.toString())
    print("Their offspring has a gene expression of:")
    print(gnm3.toString())


test()
