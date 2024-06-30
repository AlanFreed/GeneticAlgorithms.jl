module testGenes

#= First you must load the genetic algorithm.

using Pkg
Pkg.add(url = "https://github.com/AlanFreed/GeneticAlgorithms.jl")
=#

using
    GeneticAlgorithms

export
    run

function run()
    println("Genes created from boolean expressions:")
    dom = Gene(dominant)
    rec = Gene(recessive)
    s = string("A dominant  gene has an expression of ", toString(dom), ".")
    println(s)
    s = string("A recessive gene has an expression of ", toString(rec), ".")
    println(s)
    println("Genes created for gene expressions:")
    dom = Gene(Expression(dominant))
    rec = Gene(Expression(recessive))
    s = string("A dominant  gene has an expression of ", toString(dom), ".")
    println(s)
    s = string("A recessive gene has an expression of ", toString(rec), ".")
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
    println("Ten arbitrarily assigned genes are:")
    for i in 1:10
        gene = Gene()
        s = string("   ", toString(gene), ", which is ")
        if isDominant(gene)
            s = string(s, "dominant,  ")
        end
        if isRecessive(gene)
            s = string(s, "recessive, ")
        end
        mutate!(gene, 0.5)
        s = string(s, "mutates to ", toString(gene))
        println(s)
    end
    c = copy(Gene(dominant))
    s = string("A copy of a dominant expression is ", toString(c))
    println(s)
    d = deepcopy(Gene(recessive))
    s = string("while a deep copy of a recessive expression is ", toString(d))
    println(s)

    println()
    println("If these answers make sense, then this test passes.")
end # run

end # testGenes
