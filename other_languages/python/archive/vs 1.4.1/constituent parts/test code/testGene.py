# -*- coding: utf-8 -*-
#!/usr/bin/env python3
import numpy as np
from genes import Gene
"""
Created on Mon Jul 10 08:46:11 2017

@author: afreed
"""

def test():
    dom = Gene(1)
    res = Gene(0)
    gen = Gene()
    prob = 0.25
    print("A dominant  gene has an expression of {}.".format(dom.toString()))
    print("A resessive gene has an expression of {}.".format(res.toString()))
    print("An arbitrarily assigned expression of {}.".format(gen.toString()))
    if dom.isEqualTo(res):
        print("They are the same.")
    else:
        print("They are distinct gene expressions.")
    print("Setting the probability of mutation at {} one finds:".format(prob))
    g = dom.copy()
    for i in range(0,10):
        print("   {}: {}".format(i, g.toString()))
        g.mutate(prob)
    print("Popping the dominant  gene returns {}.".format(dom.pop()))
    g.push(res.pop())
    print("Pushing the resessive gene assigns {}.".format(g.pop()))
    

np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
test()