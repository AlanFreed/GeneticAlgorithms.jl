module testGenes

#= First you must load physical fields and the genetic algorithm.

using Pkg
Pkg.add(url = "https://github.com/AlanFreed/PhysicalFields.jl")
Pkg.add(url = "https://github.com/AlanFreed/GeneticAlgorithms.jl")
=#

using
    GeneticAlgorithms

import
    PhysicalFields as PF

export
    run

function run()
    println("Genes created from boolean expressions:")
    dom = Gene(dominant)
    rec = Gene(recessive)
    s = string("A dominant  gene has an expression of ", toBinaryString(dom), ".")
    println(s)
    s = string("A recessive gene has an expression of ", toBinaryString(rec), ".")
    println(s)
    println("Getting their expressions, they return:")
    s = string("   ", get(dom), "  for the dominant gene, and")
    println(s)
    s = string("   ", get(rec), " for the recessive gene.")
    println(s)
    println("Setting the latter to become the former, one gets")
    set!(rec, dominant)
    s = string("   ", get(rec))
    println(s)
    println("Ten arbitrarily assigned genes with a 50% mutation rate are:")
    for i in 1:10
        gene = Gene()
        s = string("   ", toBinaryString(gene), ", which is ")
        if isDominant(gene)
            s = string(s, "dominant,  ")
        end
        if isRecessive(gene)
            s = string(s, "recessive, ")
        end
        mutate!(gene, 0.5)
        s = string(s, "mutates to ", toBinaryString(gene))
        println(s)
    end
    c = copy(Gene(dominant))
    s = string("A copy of a dominant expression is ", toBinaryString(c), ".")
    println(s)

    # test persistence
    println()
    mydir = string(pwd(), "/files/")
    myfile = "testGenes.json"
    json_stream = PF.openJSONWriter(mydir, myfile)
    toFile(dom, json_stream)
    toFile(rec, json_stream)
    PF.closeJSONStream(json_stream)
    json_stream = PF.openJSONReader(mydir, myfile)
    gene1 = fromFile(Gene, json_stream)
    gene2 = fromFile(Gene, json_stream)
    PF.closeJSONStream(json_stream)
    if gene1 == dom && gene2 == rec
        println("A test of writing to/reading from a JSON file passed.")
    else
        println("A test of writing to/reading from a JSON file failed.")
    end
    
    println()
    println("If these answers make sense, then this test passes.")
end # run

end # testGenes
