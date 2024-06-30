module testChromosomes

#= First you must load the genetic algorithm.

using Pkg
Pkg.add(url = "https://github.com/AlanFreed/GeneticAlgorithms.jl")
=#

using
    GeneticAlgorithms

export
    run



function run()
    minX = 1.0
    maxX = 5.0
    prbM = 0.01
    prbX = 0.85
    sigF = 4
    crm1 = Chromosome(minX, maxX, sigF)
    crm2 = Chromosome(minX, maxX, sigF)
    crm3 = crossover(crm1, crm2, prbM, prbX)
    s = string("These chromosomes are made of ", crm1.genes, " genes.\n")
    s = string(s, "First  chromosome has gene expression ", toString(crm1), ".\n")
    s = string(s, "Second chromosome has gene expression ", toString(crm2), ".\n")
    s = string(s, "Their offspring has a gene expression ", toString(crm3), ".\n")
    print(s)
    if crm1 == crm3 || crm2 == crm3
        println("   The child is a clone.")
    else
        println("   The child is unique.")
    end
    println("Each of these chromosomes represents 1 of ", crm1.expressions, 
        " possible expressions.")
    println("These chromosomes decode as:")
    x1 = decode(crm1)
    x2 = decode(crm2)
    x3 = decode(crm3)
    s  = ""
    s  = string(s, "   For the first:  ", x1, "\n")
    s  = string(s, "   For the second: ", x2, "\n")
    s  = string(s, "   For the third:  ", x3, "\n")
    s  = string(s, "And then encode as:\n")
    encode!(crm1, x1)
    encode!(crm2, x2)
    encode!(crm3, x3)
    s = string(s, "   For the first:  ", toString(crm1), "\n")
    s = string(s, "   For the second: ", toString(crm2), "\n")
    s = string(s, "   For the third:  ", toString(crm3))
    println(s)

    println()
    println("If these answers make sense, then this test passes.")
end # run

end # testChromosomes
