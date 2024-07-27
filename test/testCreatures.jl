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
    minX = Vector{Real}(undef, 3)
    minX[1] = -1.0
    minX[2] =  3.141
    minX[3] =  1.0
    maxX = Vector{Real}(undef, 3)
    maxX[1] =  1.0
    maxX[2] =  3.141
    maxX[3] = 10.0
    constrained = Vector{Tuple}(undef, 0)
    optX = Vector{Real}(undef, 3)
    optX[1] = 0.1
    optX[2] = 3.141
    optX[3] = 5.0
    prbM = 0.01
    prbX = 0.85
    sigF = 4
    cre1 = procreate(minX, maxX, constrained, sigF)
    cre2 = alien(optX, minX, maxX, constrained, sigF)
    cre3 = conceive(cre1, cre2, constrained, prbM, prbX)
    print("")
    println("First creature has a gene expression of:")
    println(tostring(cre1))
    println("Second creature has a gene expression of:")
    println(tostring(cre2))
    println("Their offspring has a gene expression of:")
    println(tostring(cre3))
    if cre1 == cre3 || cre2 == cre3
        println("The child is a clone.")
    else
        println("The child is unique.")
    end
    println("Their created fitness were:")
    println("   For the procreated: ", get(cre1.fitness))
    println("   For the alien:      ", get(cre2.fitness))
    println("   For the conceived:  ", get(cre3.fitness))
    println("Their assigned fitness were:")
    set!(cre1.fitness, 1.0)
    set!(cre2.fitness, 2.0)
    set!(cre3.fitness, 3.0)
    println("   For the procreated: ", get(cre1.fitness))
    println("   For the alien:      ", get(cre2.fitness))
    println("   For the conceived:  ", get(cre3.fitness))
    println("While their parameters were:")
    println("   For the procreated: ", phenotypes(cre1))
    println("   For the alien:      ", phenotypes(cre2))
    println("   For the conceived:  ", phenotypes(cre3))
end # run

end # testCreatures

