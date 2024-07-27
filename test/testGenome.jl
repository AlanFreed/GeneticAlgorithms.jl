module testGenome

#= First you must load the genetic algorithm.

using Pkg
Pkg.add(url = "https://github.com/AlanFreed/GeneticAlgorithms.jl")
=#

using
    GeneticAlgorithms

export
    run

function run()
    minX = Vector{Real}(undef, 3)
    minX[1] =  0.0
    minX[2] = -1.0
    minX[3] = -5.0
    maxX = Vector{Real}(undef, 3)
    maxX[1] =  1.0
    maxX[2] =  1.0
    maxX[3] = -1.0
    prbM = 0.01
    prbX = 0.85
    sigF = 4
    gnm1 = Genome(minX, maxX, sigF)
    gnm2 = Genome(minX, maxX, sigF)
    gnm3 = crossover(gnm1, gnm2, prbM, prbX)
    println("These genomes have ", gnm1.genes, " genes and ", gnm1.chromosomes,
        " chromosomes.")
    println()
    println("The first genome has a gene expression of:")
    println(tostring(gnm1))
    println("The second genome has a gene expression of:")
    println(tostring(gnm2))
    println("Their offspring has a gene expression of:")
    println(tostring(gnm3))
    if gnm1 == gnm3 || gnm2 == gnm3
        println("The child is a clone.")
    else
        println("The child is unique.")
    end
    println("These genomes decode as:")
    x1 = decode(gnm1)
    x2 = decode(gnm2)
    x3 = decode(gnm3)
    println("For the first:  ", x1)
    println("For the second: ", x2)
    println("For the third:  ", x3)
    println("And then encode as vectors:")
    encode!(gnm1, x1)
    encode!(gnm2, x2)
    encode!(gnm3, x3)
    println("First genome has a gene expression of:")
    println(tostring(gnm1))
    println("Second genome has a gene expression of:")
    println(tostring(gnm2))
    println("Their offspring has a gene expression of:")
    println(tostring(gnm3))
end # run

end # testGenome

