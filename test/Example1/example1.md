# Genetic Algorithm Example Problem

This problem was given to my engineering students at SVSU for their final project in the winter semester of 2011.

The author [1] derived a hypo-elastic theory for describing the passive response of the soft tissue in lung, a.k.a. the parenchyma. In that document he solved the boundary value problem (BVP) of simple extension, i.e., the stretching of a rod, as it applies to this material model. Although this is not an experiment typically used to characterize lung, nor is it representative of natural lung response, nevertheless it is the most prevalent experiment used to characterize almost all materials.  The outcome
was a Riccati differential equation with constant coefficients whose solution is
$$
   \sigma = E \frac{\sinh ( a \epsilon )}{a \cosh ( a \epsilon ) -
   b \sinh (a \epsilon)} ,
$$
with material constants of
$$
    E > 0, \quad 0 < b < a ,
$$
where $\epsilon$ is the strain (viz., the input or control) and $\sigma$ is the stress (i.e., the output or response), while $a$, $b$, and $E$ are material parameters with $E$ being the elastic or Young's modulus.  Parameters $a$ and $b$ are dimensionless, while $E$ has units of stress (force per unit area).

The hypothetical experimental data sets supplied to my students are listed below:

$\epsilon = \{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,$

$\phantom{\epsilon = \{} 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0\}$

$\sigma = \{0.050, 0.111, 0.193, 0.290, 0.349, 0.450, 0.559, 0.622, 0.744, 0.835,$

$\phantom{\sigma = \{} 1.032, 1.144, 1.266, 1.396, 1.409, 1.494, 1.625, 1.675, 1.700, 1.710\}$

less an initial condition of $(\epsilon_0, \sigma_0) = (0,0)$ that is trivially satisfied. 

The minimum and maximum bounds on these parameters were assigned values of

$\{0.1, 0.1, 0.1\}$ and $\{1.5, 1.5, 1.5\}$

which just happen to be the same for all three parameters for this problem.

The time required to run this optimization problem was on the order of a few seconds.

# Report

## Graph

A graphical view of the fitted model is shown below.
![localImage](./GA_example1.png)
A fitting of model parameters for example 1, with the fit parameters being:

a = 1.00030, b = 0.74755 and E = 0.50882 kg/(m⋅s²).

## A Typical Report (translated to markdown)

For generation 1 of 52:

Fitness statistics for generation 1 with a population size of 116:

    optimum fitness  2.21132E+01,
    median           4.66315,
    arithmetic mean  5.73607,
    std deviation    4.26777,
    skewness         1.33068,
    excess kurtosis  1.95835.

The genome from the most fit creature:

1. 1111110011111010110000100101011111100000
2. 0101100111001100000110001011100000011011
3. 0101011010111101011101110010010101000101

Lists for [parameter_min, parameter_best, parameter_max]:

1. a ∈ [0.10000, 1.02244, 1.50000] 
2. b ∈ [0.10000, 0.70447, 1.50000] 
3. E ∈ [0.10000, 0.65145, 1.50000] kg/(m⋅s²)

Values for best parameter ± error, where a RMSE was computed wrt best values. Data were the parameters from all 116 adults living in generation 1.

1. a = 1.02244 ± 0.31562
2. b = 0.70447 ± 0.33823
3. E = 0.65145 ± 0.44040 kg/(m⋅s²)

For generation 52 of 52:

Fitness statistics for generation 52 with a population size of 116:

    optimum fitness  8.79553E+01,
    median           8.78021E+01,
    arithmetic mean  7.78090E+01, 
    std deviation    2.50535E+01,
    skewness        -2.29031,
    excess kurtosis  3.49048. 

To aid in assessing goodness of fit, a suite of correlation coefficients
between experimental response and model prediction are provided.

    Pearson's linear correlation coefficient:
      r = 0.99849,
    Spearman's monotonic correlation coefficient:
      ρ = 1.00000,
    Kendall's monotonic correlation coefficient:
      τ = 1.00000,
    and Chatterjee's nonlinear correlation coefficient:
      ξ = 0.85714.

The genome from the most fit creature:

1. 1111011011110000000100110001101101110010
2. 0100110101011100101111110101001101010010
3. 0110111110100001111100101101010011001000

Lists for [parameter_min, parameter_best, parameter_max]:

1. a ∈ [0.10000, 1.00030, 1.50000] 
2. b ∈ [0.10000, 0.74755, 1.50000] 
3. E ∈ [0.10000, 0.50882, 1.50000] kg/(m⋅s²)

Values for best parameter ± error, where a RMSE was computed wrt best values. Data were the parameters from all 116 adults living in generation 52.

1. a = 1.00030 ± 0.04458
2. b = 0.74755 ± 0.02383
3. E = 0.50882 ± 0.08596 kg/(m⋅s²)

**NOTE**: This solution space has multiple local minima, any of which the solver may converge on, so you will likely want to run the genetic algorithm several times to find the *best* one. Increasing the number of significant figures ought to help, too.  This is a probabilistic algorithm, so successive runs will most likely lead to different answers; answers that distribute about the actual solution.

# Reference

1) Freed, A.D. and Einstein, D.R., "Viscoelastic model for lung
  parenchyma for multi-scale modeling of respiratory system. Phase I:
  Hypo-elastic model for CFD implementation," Technical Report PNNL-20340,
  Pacific Northwest National Laboratory, Richland, WA, 2011.

