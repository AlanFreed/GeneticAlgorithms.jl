"""
An implementation of the parent type

    abstract type Species end

is where the user interfaces with the engine of this genetic algorithm, i.e.,
this module is what the user must program to solve their optimization problem.

To create such an implementation, one must code in the following structure:

struct mySpecies <: Species
    nExp::Int64
    nPts::Vector{Int64}
    nCtl::Vector{Int64}
    nRes::Vector{Int64}
    expC::Vector{Matrix{Float64}}
    expR::Vector{Matrix{Float64}}
    modR::Vector{Matrix{Float64}}
end # mySpecies

Parent class Species is to be inherited by all applications, e.g., by mySpecies,
where a model is to be parameterized by a genetic algorithm.  Here is where the
user creates their 'model' whose 'parameters' a user seeks to quantify through
this optimization scheme.

The model is to assign its parameter array so that it indexes from 2.  This is
because index 1 is reserved for internal use by the optimizer to determine its
optimal objective function parameter.  This is a vital feature of this scheme.

Constructor

    s = mySpecies(expTypes, ptsInExp, ctrlsInExp, respsInExp)

where

    expTypes    number of experiment types with separate data sets
    ptsInExp    number of datum points in experiment[index]
    ctlInExp    number of controlled variables in experiment[index]
    rspInExp    number of  response  variables in experiment[index]

with ptsInExp, ctlInExp, rspInExp being integer arrays.

Attributes (fields of the type) accessed by the genetic algorithm are:

    nExp    integer:         number of different experiment types
    nPts    integer array:   number of data points in each experiment type
    nCtl    integer array:   number of control variables per experiment type
    nRsp    integer array:   number of response variables per experiment type
    expC    real array:      vector of matrices of experimental control values
                             sized at exp[i] X [nCtl[i] X nPts[i]]
                             with i indexing from 1 to nExp
    expR    real array:      vector of matrices of experimental response values
                             sized at exp[i] X [nRsp[i] X nPts[i]]
                             with i indexing from 1 to nExp
    stdE    real array:      vector of vectors for standard deviations of expR
                             sized at exp[i] X nRsp[i]
                             with i indexing from 1 to nExp
    modR    real array:      vector of matrices of model response values
                             sized at exp[i] X [nRsp[i] X nPts[i]]
                             with i indexing from 1 to nExp

Methods

The following methods are called prior to calling the genetic algorithm.

    data = newExpCtlData(s)
        data                creates a new array of matrices to be filled in by
                            a user, e.g., data[i][j,k], where:
                            i)  vector indexes over the number of experiments
                            ii) matrix has rows that index over the number
                                of controlled variables for the specified
                                experiment, and columns that index over the
                                number of data points for that experiment
    assignExpCtlData!(s, data)
        data                assigns a filled-in data structure to field 'expC'
    data = newExpRspData(s)
        data                creates a new array of matrices to be filled in by
                            a user, e.g., data[i][j,k], where:
                            i)  vector indexes over the number of experiments
                            ii) matrix has rows that index over the number
                                of response variables for the specified
                                experiment, and columns that index over the
                                number of data points for that experiment
    assignExpRspData!(s, data)
        data                assigns a filled-in data structure to field 'expR'

The following method is called from within the genetic algorithm.

    runModel(s, parameters) solves the user's model using the parameters given
                            and places the response results within field 'modR'
        parameters          array of model parameters managed by the genetic
                            algorithm: start indexing your model's parameters
                            with index 2, as index 1 is reserved for an internal
                            parameter within the objective function itself
"""
abstract type Species end
