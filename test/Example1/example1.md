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

# Reference

1) Freed, A.D. and Einstein, D.R., "Viscoelastic model for lung
  parenchyma for multi-scale modeling of respiratory system. Phase I:
  Hypo-elastic model for CFD implementation," Technical Report PNNL-20340,
  Pacific Northwest National Laboratory, Richland, WA, 2011.

