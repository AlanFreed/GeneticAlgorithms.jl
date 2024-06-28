


This module originated from courses in numerical methods that the author (Alan Freed) taught at both Saginaw Valley State University (SVSU) and at Texas A&M University (TAMU).

To use this module, you will need to add the following repository to your project:

```
using Pkg
Pkg.add(url = "https://github.com/AlanFreed/GeneticAlgorithms.jl")
```

# GeneticAlgorithms

Details of the genetic algorithm are buried within the library.  They do not permeate up to the interface.  The author drew heavily on the books of Goldberg cited below in his construction of this genetic algorithm.

At the core of this genetic algorithm are haploid genes that can admit two expressions, passive and recessive, represented by a binary value.  At a very low probability of occurrence, a gene may mutate from passive to active or vice versa.

A chromosome is an array of genes assembled as a Grey binary number.  Grey numbers are selected over regular binary numbers so that mutation events have a more refined effect.  Each chromosome represents a model parameter whose optimal value is being sought.  A decoder/encoder pair maps a floating point real into a Grey binary array and back again.  These are one-to-one mappings.  

A genome is a collection of chromosomes, one for each parameter being sought.  It contains the genetic information.

A creature is comprised of a genome and a fitness that associates with it.  Fitness is taken to be a range-weighted harmonic sum of the root mean squared errors between model predictions acquired over a collection of experiments.  

A colony is a population of creatures, usually numbering in the hundreds.  To evolve a colony from one generation into the next, creatures from the existing generation are picked at random for mate selection via tournament play.  From two tournaments, a pair of parent creatures get selected for mating with the outcome being a child creature that will move onto the next generation.  Only the elite creature from an existing population moves onto the next generation.  All other creatures perish at the moment when a new population of children from the current generation become the adults of the next generation.  Mating involves a crossover event with a high probability of occurrence between like pairs of chromosomes from two parents.  Where a pair of chromosomes are to be split and their genetic information swapped (i.e., a crossover event) which is a probabilistic outcome.  The end result is that a child's chromosome strands are comprised of gene segments from both parents.  

The original generation is procreated.  Every gene is assigned randomly with 50-50 odds in that generation.  This provides a dispersed population over the window of admissible parametric values, a Monte Carlo sampling.  To aid in maintaining a diverse population, a procreated immigrant or two may get introduced into a population at each generation at a low probability of occurrence.  Mutation also aids in maintaining diversity.  No clones are allowed in this implementation of a genetic algorithm.  Instead of evolving a population of creatures to an existence where all creatures are clones of one another, this implementation of a genetic algorithm determines its convergence criteria and population size according to algorithms derived by Goldberg (2002) from probability theory.

# References

	1. Goldberg, D.E., Genetic Algorithms in Search, Optimization, and Machine Learning, Addison-Wesley, Boston, 1989.
    2. Goldberg, D.E., The Design of Innovation: Lessons learned from and for competent genetic algorithms.  In: Genetic algorithms and evolutionary computation, Vol. 7, Klewer, Boston, 2002.    
 
# Version History

# Version 0.1.0

The original version, dated June 2024, is a port from the author's prior codes written in Zonnon, Pascal and Python for classes that he taught in numerical methods to mechanical engineers when he was a professor at SVSU (2007-2014) and TAMU (2014-2021).
