module testChromosomes

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
    minX = PF.PhysicalScalar(PF.CENTIMETER)
    PF.set!(minX, 1.0)
    maxX = PF.PhysicalScalar(PF.CENTIMETER)
    PF.set!(maxX, 5.0)
    prbM = 0.01
    prbX = 0.85
    sigF = 5
    crm1 = Chromosome(minX, maxX, sigF)
    crm2 = Chromosome(minX, maxX, sigF)
    crm3 = crossover(crm1, crm2, prbM, prbX)
    s = string("These chromosomes are made of ", crm1.genes, " genes.\n")
    s = string(s, "First  chromosome has gene expression ", toBinaryString(crm1), ".\n")
    s = string(s, "Second chromosome has gene expression ", toBinaryString(crm2), ".\n")
    s = string(s, "Their offspring has a gene expression ", toBinaryString(crm3), ".\n")
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
    s  = string(s, "   For the first:  ", PF.toString(x1), "\n")
    s  = string(s, "   For the second: ", PF.toString(x2), "\n")
    s  = string(s, "   For the third:  ", PF.toString(x3), "\n")
    s  = string(s, "And then encode as:\n")
    encode!(crm1, x1)
    encode!(crm2, x2)
    encode!(crm3, x3)
    s = string(s, "   For the first:  ", toBinaryString(crm1), "\n")
    s = string(s, "   For the second: ", toBinaryString(crm2), "\n")
    s = string(s, "   For the third:  ", toBinaryString(crm3))
    println(s)

    # test persistence
    println()
    mydir = string(pwd(), "/files/")
    myfile = "testChromosomes.json"
    json_stream = PF.openJSONWriter(mydir, myfile)
    toFile(crm1, json_stream)
    toFile(crm2, json_stream)
    toFile(crm3, json_stream)
    PF.closeJSONStream(json_stream)
    json_stream = PF.openJSONReader(mydir, myfile)
    x1 = fromFile(Chromosome, json_stream)
    x2 = fromFile(Chromosome, json_stream)
    x3 = fromFile(Chromosome, json_stream)
    PF.closeJSONStream(json_stream)
    if crm1 == x1 && crm2 == x2 && crm3 == x3
        println("A test of writing to/reading from a JSON file passed.")
    else
        println("A test of writing to/reading from a JSON file failed.")
    end
    println()
    println("If these answers make sense, then this test passes.")
end # run

end # testChromosomes
