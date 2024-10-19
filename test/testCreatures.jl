module testCreatures

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
    minX[1] = PF.PhysicalScalar(-1.0, PF.CGS_STRESS)
    minX[2] = PF.PhysicalScalar(3.141, PF.CGS_STRAIN_RATE)
    minX[3] = PF.PhysicalScalar(1.0, PF.CGS_STRESS)
    maxX = Vector{PF.PhysicalScalar}(undef, 3)
    maxX[1] = PF.PhysicalScalar(1.0, PF.CGS_STRESS)
    maxX[2] = PF.PhysicalScalar(3.141, PF.CGS_STRAIN_RATE)
    maxX[3] = PF.PhysicalScalar(10.0, PF.CGS_STRESS)
    optX = Vector{PF.PhysicalScalar}(undef, 3)
    optX[1] = PF.PhysicalScalar(0.1, PF.CGS_STRESS)
    optX[2] = PF.PhysicalScalar(3.141, PF.CGS_STRAIN_RATE)
    optX[3] = PF.PhysicalScalar(5.0, PF.CGS_STRESS)
    cnst = Vector{Tuple{Int,Int}}(undef, 0)
    prbM = 0.01
    prbX = 0.95
    sigF = 4
    cre1 = procreate(minX, maxX, cnst, sigF)
    cre2 = alien(optX, minX, maxX, cnst, sigF)
    cre3 = conceive(cre1, cre2, cnst, prbM, prbX)
    print("")
    println("First creature has a gene expression of:")
    println(toBinaryString(cre1))
    println("Second creature has a gene expression of:")
    println(toBinaryString(cre2))
    println("Their offspring has a gene expression of:")
    println(toBinaryString(cre3))
    if cre1 == cre3 || cre2 == cre3
        println("The child is a clone.")
    else
        println("The child is unique.")
    end
    println()
    println("Their created fitness were:")
    println("For the procreated:")
    println("   ", PF.toString(cre1.fitness))
    println("For the alien:")
    println("   ", PF.toString(cre2.fitness))
    println("For the conceived:")
    println("   ", PF.toString(cre3.fitness))
    println()
    println("Their assigned fitness were:")
    PF.set!(cre1.fitness, 1)
    PF.set!(cre2.fitness, 2)
    PF.set!(cre3.fitness, 3)
    println("For the procreated:")
    println("   ", PF.toString(cre1.fitness))
    println("For the alien:")
    println("   ", PF.toString(cre2.fitness))
    println("For the conceived:")
    println("   ", PF.toString(cre3.fitness))
    println()
    println("While their parameters were:")
    println("For the procreated: ")
    println(toString(cre1.DNA))
    println("For the alien:")
    println(toString(cre2.DNA))
    println("For the conceived:")
    println(toString(cre3.DNA))
    # test persistence
    println()
    mydir = string(pwd(), "/files/")
    myfile = "testCreatures.json"
    json_stream = PF.openJSONWriter(mydir, myfile)
    toFile(cre1, json_stream)
    toFile(cre2, json_stream)
    toFile(cre3, json_stream)
    PF.closeJSONStream(json_stream)
    json_stream = PF.openJSONReader(mydir, myfile)
    red1 = fromFile(Creature, json_stream)
    red2 = fromFile(Creature, json_stream)
    red3 = fromFile(Creature, json_stream)
    PF.closeJSONStream(json_stream)
    if cre1 == red1 && cre2 == red2 && cre3 == red3
        println("A test of writing to/reading from a JSON file passed.")
    else
        println("A test of writing to/reading from a JSON file failed.")
    end
    println()
    println("If these answers make sense, then this test passes.")
end # run

end # testCreatures

