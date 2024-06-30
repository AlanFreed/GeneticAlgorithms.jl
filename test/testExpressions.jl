module testExpressions

#= First you must load the genetic algorithm.

using Pkg
Pkg.add(url = "https://github.com/AlanFreed/GeneticAlgorithms.jl")
=#

using
    GeneticAlgorithms

export
    run

function run()
    dom = Expression(dominant)
    rec = Expression(recessive)
    s = string("A dominant  gene has an expression of ", toString(dom), ".")
    println(s)
    s = string("A recessive gene has an expression of ", toString(rec), ".")
    println(s)
    println("Getting their expressions, they return:")
    s = string("   ", get(dom), "  for the dominant expression, and")
    println(s)
    s = string("   ", get(rec), " for the recessive expression.")
    println(s)
    println("Setting the latter to become the former, one gets")
    set!(rec, dominant)
    s = string("   ", get(rec))
    println(s)
    println("Ten arbitrarily assigned gene expressions are:")
    for i in 1:10
        exp = Expression()
        s = string("   ", toString(exp), " which is ")
        if isDominant(exp)
            s = string(s, "dominant")
        end
        if isRecessive(exp)
            s = string(s, "recessive")
        end
        println(s)
    end
    c = copy(Expression(dominant))
    s = string("A copy of a dominant expression is ", toString(c))
    println(s)
    d = deepcopy(Expression(recessive))
    s = string("while a deep copy of a recessive expression is ", toString(d))
    println(s)

    println()
    println("If these answers make sense, then this test passes.")
end # run

end # testExpressions
