# Genetic Algorithms

This module originated from courses in numerical methods that the author (Alan Freed) taught at both Saginaw Valley State University (SVSU) and at Texas A&M University (TAMU).

To use this module, you will need to add the following repository to your project:

```
using Pkg
Pkg.add(url = "https://github.com/AlanFreed/GeneticAlgorithms.jl")
```

# Overview

Goldberg [1] tells his readers what genetic algorithms (GAs) are, conceptually, in the opening paragraph of his book:

> Genetic algorithms are search algorithms based on the mechanics of natural selection and natural genetics. They combine survival of the fittest among string structures with a structured yet randomized information exchange to form a search algorithm with some of the innovative flair of human search.  In every generation, a new set of artificial creatures (strings) is created using bits and pieces of the fittest of the old; an occasional new part is tried for good measure. While randomized, genetic algorithms are no simple random walk. They efficiently exploit historical information to speculate on new search points with expected improved performance.

GAs have in common a population of individuals, a means to determine the
fitness for each individual within the population, the pairing of individuals
for reproduction according to their fitness, the cross-fertilization of
genetic material from parents to produce offspring, and a random chance of a
mutation occurring in an off-spring's genetic code.  A new generation
replaces the existing one, and the reproduction cycle starts all over again.
This process is repeated to convergence.

This implementation of a genetic algorithm draws
heavily on the algorithms and Pascal code of Goldberg [1,2].  His books focus on what he calls a
simple genetic algorithm (SGA), which is implemented here. 
Schmitt [3] has shown that GAs can be derived
from the principle of maximum entropy governing a Markovian process.

GAs are but one of many optimization techniques that exist.  There is a
theorem in the optimization literature called the "no free lunch"
theorem.  It states that the overall performance of any optimization
algorithm, when evaluated over the set of all possible optimization
problems, is no different than any other optimization technique.  An
apparent advantage of one
algorithm over another resides strictly with its application and, quite
often, the personal taste of the user.  It has been the author's
experience that GAs work very well for parameter estimations, particularly for
non-linear models.

# Statistics

The objective functions employed herein are based upon sample statistics pertaining to experimental responses, denoted as **x** = {*x*₁, *x*₂, …, *x*ₙ}, and their corresponding predicted model responses, denoted as **ξ** = {*ξ*₁, *ξ*₂, …, *ξ*ₙ}, where *n* signifies the number of vector entries or datum points obtained from the experiment. These data are normalized by defining **X** = {*X*₁, *X*₂, …, *X*ₙ} = {*x*₁, *x*₂, …, *x*ₙ} / *max*(|**x**|) so that *X*ᵢ ∈ [-1, 1] ∀ *i*  ∈ [1, *n*]. The model responses are also normalized; specifically, **Ξ** = {*Ξ*₁, *Ξ*₂, …, *Ξ*ₙ} = {*ξ*₁, *ξ*₂, …, *ξ*ₙ} / *max*(|**x**|). (**NOTE:** The normalized model response **Ξ** normalizes **ξ** with the experimental value *max*(|**x**|), not with the model's value *max*(|**ξ**|)).

An error can therefore be defined as 

**ϵ** = **X** - **Ξ**.

## Sample statistics

The expectation *E* of a data set **X** is the first moment of its data, viz., 

*E*(**X**) = (1/N) ∑ᵢ₌₁ⁿ *X*ᵢ

The sample variance of this data set is

*VAR*(**X**) = *E*((**X** - *E*(**X**))²) = *E*(**X**²) - (*E*(**X**))² where *E*(**X**²) is the second moment (1/N) ∑ᵢ₌₁ⁿ *X*ᵢ²

while the sample covariance between two data sets is

*COV*(**X**, **Ξ**) = *E*((**X** - *E*(**X**))(**Ξ** - *E*(**Ξ**))) = *E*(**XΞ**) - *E*(**X**) *E*(**Ξ**) where *E*(**XΞ**) is a mixed moment (1/N) ∑ᵢ₌₁ⁿ *X*ᵢ*Ξ*ᵢ

which reduces to the sample variance whenever **Ξ** = **X**.

## Objective functions

A collection of objective functions that have statisical warrent are considered.

1. Minimize expectation for the magnitude of error.
    1. ϕ₁ = 1 / *E*(|**ϵ**|) = 1 / *E*(|**X** - **Ξ**|)
2. Minimize expectation for the squared error.
    1. ϕ₂ = 1 / *E*(**ϵ**²) = 1 / (*E*((**X** - **Ξ**)²) = 1 / (*VAR*(**X**) + *VAR*(**Ξ**) - 2*COV*(**X**, **Ξ**) + (*E*(**X**))² + (*E*(**Ξ**))² - 2*E*(**X**) *E*(**Ξ**))
3. Minimize variance for the error.
    1. ϕ₃ = 1 / *VAR*(**ϵ**) = 1 / (*VAR*(**X**) + *VAR*(**Ξ**) - 2*COV*(**X**, **Ξ**) )
4. Maximize covariance between experiment **X** and model **Ξ**.
    1. ϕ₄ = *COV*(**X**, **Ξ**)

Genetic algorithms seek to maximize a quality parameter referred to as *fitness*; hence, reciprocal values are used to describe the fitness of an objective function that minimizes.

# A Genetic Algorithm

At the core of this genetic algorithm are haploid genes that can admit two expressions, dominate and recessive. They are represented with a binary value.  At a very low probability of occurrence, a gene may mutate from dominate to recessive, or vice versa.

A chromosome is an array of genes assembled as a Gray binary number.  Gray numbers are selected over regular binary numbers so that mutation events have a more refined effect.  Each chromosome represents a model parameter whose optimal value is being sought.  A decoder/encoder pair maps a floating point real into a Gray binary array and back again.  These are one-to-one mappings.  

A genome is a collection of chromosomes, one for each model parameter being sought.  It contains the genetic information of a creature.

A creature is comprised of a genome and a fitness that associates with it.  Fitness is taken to be a root mean squared error between model predictions
and experimental outcomes, acquired over a collection of experiments. These errors are normalized by the standard deviation in the experimental response, which allows data between different types of experiments to be 
used for parameterization in a meaningful way.

A colony is a population of creatures, usually numbering in the hundreds.  To evolve a colony from one generation into the next, creatures from the existing generation are picked at random for mate selection via tournament play.  From two tournaments, a pair of parent creatures get selected for mating with the outcome being a child creature that will move onto the next generation.  Only the elite creature from an existing population moves onto the next generation.  All other creatures perish at the moment when a new population of children from the current generation become the adults of the next generation.  Mating involves a crossover event with a high probability of occurrence between like pairs of chromosomes from two parents.  During such an event, a pair of chromosomes are split and their genetic information swapped (i.e., a crossover event occurs) which is a probabilistic outcome.  The end result is that a child's chromosome strands are comprised of gene segments from both parents.  

The original generation is procreated.  Every gene is assigned an expression randomly with 50-50 odds in that generation.  This provides a dispersed population over the window of admissible parametric values, i.e., a Monte Carlo sampling.  To aid in maintaining a diverse population, a procreated immigrant or two may get introduced into a population at each generation at a low probability of occurrence.  Mutation also aids in maintaining diversity.  No clones are allowed in this implementation of a genetic algorithm.  Instead of evolving a population of creatures to an existence where all creatures are clones of one another, which is often done in genetic algorithms, this implementation of a genetic algorithm determines its convergence criteria and population size according to algorithms derived by Goldberg [2] from probability theory.

It bears repeating a phrase that commonly appears: "Correlation
does not imply causation."  In other words, good correlation must not be used to infer or suggest the existence of a causal relationship between a
model and data, i.e., it does not prove that the model is 'correct'.

## Mutable types

### Expressions

### Counters

### Variables

## Genetic types

### Genes

### Chromosomes

### Genome

### Creatures

## Algorithmic types

### Colony

### User interface

### Genetic algorithm

# Examples

# References

	1. Goldberg, D.E., Genetic Algorithms in Search, Optimization, and Machine Learning, Addison-Wesley, Boston, 1989.
    2. Goldberg, D.E., The Design of Innovation: Lessons learned from and for competent genetic algorithms, 2002.  In: Genetic algorithms and evolutionary computation, Vol. 7, Klewer, Boston, 2002.
    3. Schmitt, L.M., "Theory of genetic algorithms", Theoretical Computer Science, Vol. 259 (2001) pp. 1-61.
    4. Schmitt, L.M., "Theory of genetic algorithms II: Models for genetic operators over the string-tensor representation of populations and convergence to global optima for arbitrary fitness function under scaling", Theoretical Computer Science, Vol. 310 (2004) pp. 181-231.
    5. Sivanandam, S.N. and Deepa, S.N., Introduction to Genetic Algorithms, Springer, Berlin, 2008.
 
# Version History

## Version 0.1.3

Four fitness types, i.e., objective functions, have been introduced:

    1. Minimize expectation for the magnitude of error.
    2. Minimize expectation for the squared error.
    3. Minimize sample variance in the error.
    4. Maximize sample covariance between experiment and model.
    

## Version 0.1.2

A complete rewrite focusing on the use of multiple dispatch in liu of object-oriented coding.

## Version 0.1.0

The original version, dated June 2024, is a port from the author's prior codes written in Zonnon, Pascal and Python for classes that he taught in numerical methods to mechanical engineers when he was a professor at SVSU (2007-2014) and TAMU (2014-2021).
