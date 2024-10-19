module testGenome

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

function run()
    minX = Vector{PF.PhysicalScalar}(undef, 3)
    minX[1] = PF.PhysicalScalar(-5.0, PF.CGS_STRESS)
    minX[2] = PF.PhysicalScalar(-1.0, PF.CGS_STRAIN_RATE)
    minX[3] = PF.PhysicalScalar(-1.0, PF.CGS_STRESS)
    maxX = Vector{PF.PhysicalScalar}(undef, 3)
    maxX[1] = PF.PhysicalScalar(1.0, PF.CGS_STRESS)
    maxX[2] = PF.PhysicalScalar(1.0, PF.CGS_STRAIN_RATE)
    maxX[3] = PF.PhysicalScalar(1.0, PF.CGS_STRESS)
    cnst = Vector{Tuple{Int,Int}}(undef, 1)
    cnst[1] = (1, 3)
    prbM = 0.01
    prbX = 0.95
    sigF = 5
    gnm1 = Genome(minX, maxX, cnst, sigF)
    gnm2 = Genome(minX, maxX, cnst, sigF)
    gnm3 = crossover(gnm1, gnm2, prbM, prbX)
    println("These genomes have ", gnm1.genes, " genes and ", gnm1.chromosomes,
        " chromosomes.")
    println()
    println("The first genome has a gene expression of:")
    println(toBinaryString(gnm1))
    println("The second genome has a gene expression of:")
    println(toBinaryString(gnm2))
    println("Their offspring has a gene expression of:")
    println(toBinaryString(gnm3))
    if gnm1 == gnm3 || gnm2 == gnm3
        println("The child is a clone.")
    else
        println("The child is unique.")
    end
    println()
    println("These genomes decode as:")
    println()
    x1 = decode(gnm1)
    x2 = decode(gnm2)
    x3 = decode(gnm3)
    println("For the first:")
    println(toString(gnm1))
    println("For the second:")
    println(toString(gnm2))
    println("For the third:")
    println(toString(gnm3))
    println()
    println("And then encode as vectors:")
    println()
    encode!(gnm1, x1)
    encode!(gnm2, x2)
    encode!(gnm3, x3)
    println("First genome has a gene expression of:")
    println(toBinaryString(gnm1))
    println("Second genome has a gene expression of:")
    println(toBinaryString(gnm2))
    println("Their offspring has a gene expression of:")
    println(toBinaryString(gnm3))
    # test persistence
    println()
    mydir = string(pwd(), "/files/")
    myfile = "testGenomes.json"
    json_stream = PF.openJSONWriter(mydir, myfile)
    toFile(gnm1, json_stream)
    toFile(gnm2, json_stream)
    toFile(gnm3, json_stream)
    PF.closeJSONStream(json_stream)
    json_stream = PF.openJSONReader(mydir, myfile)
    x1 = fromFile(Genome, json_stream)
    x2 = fromFile(Genome, json_stream)
    x3 = fromFile(Genome, json_stream)
    PF.closeJSONStream(json_stream)
    if gnm1 == x1 && gnm2 == x2 && gnm3 == x3
        println("A test of writing to/reading from a JSON file passed.")
    else
        println("A test of writing to/reading from a JSON file failed.")
    end
    println()
    println("If these answers make sense, then this test passes.")
end # run

end # testGenome

