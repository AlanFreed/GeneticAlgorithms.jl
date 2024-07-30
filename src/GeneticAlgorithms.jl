#=
Created on Thr 14 Jun 2024
Updated on Thr 25 Jul 2024
Translated from python (code dated 09/08/2017) with enhancements for julia.
=#

"""
This julia module is a translation of a python module written for the author's
course on numerical methods at TAMU, which was a translation of a pascal module,
which was a translation of the author's original genetic algorithm written in
zonnon that the author first used to illustrate to his students at SVSU the art
of writing object-oriented code, for which genetic algorithms are well suited.

From a biologic interpretation, the genetic algorithm implemented here is a
colony of creatures that advances their quality of life (fitness) from one
generation to the next.  Class 'GeneticAlgorithm,' defined at the end of this
file, is a colony or collection of creatures whose population is sustained from
one generation to the next.  Mating between creatures occurs through a process
known as tournament play, where the most fit contestant from a random selection
of contestants is chosen for mating.  Typically, each successive generation is
more fit than its predecessor, i.e., the colony's quality improves with time.

To install this package, download the following package from its URL address:

using Pkg
Pkg.add(url = "https://github.com/AlanFreed/GeneticAlgorithms.jl")

The data structure for genetic algorithm has an interface of

struct GeneticAlgorithm
    c::Colony
end

with constructor

    ga = GeneticAlgorithm(colony)

and procedure

    run(ga, verbose)

    Function 'run' runs a solver for genetic algorithm 'ga' whose population
    size and number of generations to advance through are determined internally.
    A report is written to file for the user to read.  If 'verbose' is true,
    the default, then a page in this report is written for each generation of
    the colony; otherwise, the report only contains information pertaining to
    the final generation.
"""
module GeneticAlgorithms

using
    Statistics,
    StatsBase

import
    Printf: @sprintf

export
    # abstract type
    AbstractParameters,

    # genetic types
    Expression,
    Counter,
    Variable,
    Gene,
    Chromosome,
    Genome,
    Creature,
    Colony,

    # external constructors for type Creature
    procreate,
    alien,
    conceive,

    # algorithm types
    ExperimentalData,
    Model,
    GeneticAlgorithm,

    # operators:
    # logical
    ==, ≠, <, ≤, >, ≥,
    # arithmetic unary
    +, -,
    # arithmetic binary
    +, -, *, ÷, %, /, ^,

    # methods
    copy,
    get,
    getindex,
    set!,
    setindex!,
    tostring,

    isdominant,
    isrecessive,
    mutate!,
    decode,
    encode!,
    crossover,
    phenotypes,
    solve,
    advance_to_next_generation!,
    report,
    run,

    # constants
    dominant,
    recessive

include("Expressions.jl")

include("Counters.jl")

include("Variables.jl")

include("Genes.jl")

include("Chromosomes.jl")

include("Genomes.jl")

include("Creatures.jl")

include("UserInterface.jl")

include("Colonies.jl")

# The genetic algorithm.

struct GeneticAlgorithm
    colony::Colony

    # constructor

    function GeneticAlgorithm(colony::Colony)
        new(colony)
    end
end # GeneticAlgorithm

# Method

function run(ga::GeneticAlgorithm, verbose::Bool=true)
    println("Each ⋅ represents one generation advanced out of ", ga.colony.generations_to_convergence, " generations total.")
    print("    ")

    # Open an IO-stream to write to.
    mydir = string(pwd(), "/files")
    if !isdir(mydir)
        mkdir(mydir)
    end
    myfile = string(mydir, "/ga_report.txt")
    if isfile(myfile)
        mystream = open(myfile; lock=true, read=false, write=true, create=false, truncate=true, append=true)
    else
        mystream = open(myfile; lock=true, read=false, write=true, create=true, truncate=true, append=true)
    end
    seekstart(mystream)

    # The first generation.
    s = string("For generation ", tostring(ga.colony.generation), " of ", ga.colony.generations_to_convergence, ":\n\n")
    s = string(s, report(ga.colony))
    write(mystream, s)
    flush(mystream)

    # Run the genetic algorithm.
    for i in 2:ga.colony.generations_to_convergence-1
        print("⋅")
        advance_to_next_generation!(ga.colony)
        if verbose
            s = string("For generation ", tostring(ga.colony.generation), " of ", ga.colony.generations_to_convergence, ":\n\n")
            s = string(s, report(ga.colony))
            write(mystream, s)
            flush(mystream)
        end
    end
    println("⋅")
    advance_to_next_generation!(ga.colony)
    println("The genetic algorithm has finished.\n\n")
    s = string("For generation ", tostring(ga.colony.generation), " of ", ga.colony.generations_to_convergence, ":\n\n")
    s = string(s, report(ga.colony))
    write(mystream, s)
    flush(mystream)
    close(mystream)
    println("See file ", myfile, " for a report.")

    return nothing
end # run

end # module GeneticAlgorithms
