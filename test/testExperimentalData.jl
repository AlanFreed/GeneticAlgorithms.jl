module testExperimentalData

#= First you must load physical fields and the genetic algorithm.

using Pkg
Pkg.add(url = "https://github.com/AlanFreed/PhysicalFields.jl")
Pkg.add(url = "https://github.com/AlanFreed/GeneticAlgorithms.jl")
=#

import
    PhysicalFields as PF

using
    GeneticAlgorithms

export
    run

function myData()::GeneticAlgorithms.ExperimentalData
    n_exp = 2
    pairs = 2
    n_dat = 3
    
    experiments     = n_exp
    conjugate_pairs = Vector{Int}(undef, n_exp)
    data_points     = Vector{Int}(undef, n_exp)
    for exp in 1:n_exp
        conjugate_pairs[exp] = pairs
        data_points[exp]     = n_dat
    end
    
    control  = Vector{Vector{PF.ArrayOfPhysicalScalars}}(undef, n_exp)
    response = Vector{Vector{PF.ArrayOfPhysicalScalars}}(undef, n_exp)
    for exp in 1:n_exp
        control[exp]  = Vector{PF.ArrayOfPhysicalScalars}(undef, pairs)
        response[exp] = Vector{PF.ArrayOfPhysicalScalars}(undef, pairs)
        for pair in 1:pairs
            control[exp][pair]  = PF.ArrayOfPhysicalScalars(n_dat, PF.STRAIN)
            response[exp][pair] = PF.ArrayOfPhysicalScalars(n_dat, PF.STRESS)
        end
    end
    control[1][1][1] = PF.PhysicalScalar(0.1, PF.STRAIN)
    control[1][1][2] = PF.PhysicalScalar(0.2, PF.STRAIN)
    control[1][1][3] = PF.PhysicalScalar(0.3, PF.STRAIN)
    control[1][2][1] = PF.PhysicalScalar(0.4, PF.STRAIN)
    control[1][2][2] = PF.PhysicalScalar(0.5, PF.STRAIN)
    control[1][2][3] = PF.PhysicalScalar(0.6, PF.STRAIN)
    control[2][1][1] = PF.PhysicalScalar(0.7, PF.STRAIN)
    control[2][1][2] = PF.PhysicalScalar(0.8, PF.STRAIN)
    control[2][1][3] = PF.PhysicalScalar(0.9, PF.STRAIN)
    control[2][2][1] = PF.PhysicalScalar(1.0, PF.STRAIN)
    control[2][2][2] = PF.PhysicalScalar(1.1, PF.STRAIN)
    control[2][2][3] = PF.PhysicalScalar(1.2, PF.STRAIN)
    
    response[1][1][1] = PF.PhysicalScalar(0.050, PF.STRESS)
    response[1][1][2] = PF.PhysicalScalar(0.111, PF.STRESS)
    response[1][1][3] = PF.PhysicalScalar(0.193, PF.STRESS)
    response[1][2][1] = PF.PhysicalScalar(0.290, PF.STRESS)
    response[1][2][2] = PF.PhysicalScalar(0.349, PF.STRESS)
    response[1][2][3] = PF.PhysicalScalar(0.450, PF.STRESS)
    response[2][1][1] = PF.PhysicalScalar(0.559, PF.STRESS)
    response[2][1][2] = PF.PhysicalScalar(0.622, PF.STRESS)
    response[2][1][3] = PF.PhysicalScalar(0.744, PF.STRESS)
    response[2][2][1] = PF.PhysicalScalar(0.835, PF.STRESS)
    response[2][2][2] = PF.PhysicalScalar(1.032, PF.STRESS)
    response[2][2][3] = PF.PhysicalScalar(1.144, PF.STRESS)
    
    return ExperimentalData(experiments, conjugate_pairs,
                            data_points, control, response)
end # myData

function run()
    mydata = myData()
    copydata = copy(mydata)
    if copydata == mydata
        println("Method copy returned an accurate copy.")
    else
        println("Method copy did *not* return an accurate copy.")
    end
    # test persistence
    println()
    mydir = string(pwd(), "/files/")
    myfile = "testExpData.json"
    json_stream = PF.openJSONWriter(mydir, myfile)
    toFile(mydata, json_stream)
    PF.closeJSONStream(json_stream)
    json_stream = PF.openJSONReader(mydir, myfile)
    data1 = fromFile(ExperimentalData, json_stream)
    PF.closeJSONStream(json_stream)
    if mydata == data1
        println("A test of writing to/reading from a JSON file passed.")
    else
        println("A test of writing to/reading from a JSON file failed.")
    end
    
    println()
    println("If these answers make sense, then this test passes.")
end # run

end # testExperimentalData
