module testCreatures

#= First you must load the genetic algorithm.

using Pkg
Pkg.add(url = "https://github.com/AlanFreed/GeneticAlgorithms.jl")
=#

using
    GeneticAlgorithms

export
    run

function run()
    minX = [-1.0, 3.141,  1.0]
    optX = [ 0.1, 3.141,  5.0]
    maxX = [ 1.0, 3.141, 10.0]
    prbM = 0.01
    prbX = 0.85
    sigF = 4
    cre1 = procreate(minX, maxX, sigF)
    cre2 = alien(optX, minX, maxX, sigF)
    cre3 = conceive(cre1, cre2, prbM, prbX)
    print("")
    println("First creature has a gene expression of:")
    println(toString(cre1))
    println("Second creature has a gene expression of:")
    println(toString(cre2))
    println("Their offspring has a gene expression of:")
    println(toString(cre3))
    if cre1 == cre3 || cre2 == cre3
        println("The child is a clone.")
    else
        println("The child is unique.")
    end
    println("Their created fitness were:")
    println("   For the procreated: ", getFitness(cre1))
    println("   For the alien:      ", getFitness(cre2))
    println("   For the conceived:  ", getFitness(cre3))
    println("Their assigned fitness were:")
    setFitness!(cre1, 1.0)
    setFitness!(cre2, 2.0)
    setFitness!(cre3, 3.0)
    println("   For the procreated: ", getFitness(cre1))
    println("   For the alien:      ", getFitness(cre2))
    println("   For the conceived:  ", getFitness(cre3))
    println("While their parameters were:")
    println("   For the procreated: ", getParameters(cre1))
    println("   For the alien:      ", getParameters(cre2))
    println("   For the conceived:  ", getParameters(cre3))
end # run

end # testCreatures

