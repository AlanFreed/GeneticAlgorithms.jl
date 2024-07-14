#=
Created on Thr 14 Jun 2024
Updated on Sun 14 Jul 2024
Translated from Python (code dated 09/08/2017) and enhanced for Julia.
=#

"""
This module was taken from the author's course on numerical methods at TAMU.

From a biologic interpretation, the genetic algorithm implemented here is a
colony of creatures that advances their quality of life (fitness) from one
generation to the next.  Class 'GeneticAlgorithm,' defined at the end of this
file, is a colony or collection of creatures whose population is sustained from
one generation to the next.  Mating between creatures occurs through a process
known as tournament play, where the most fit contestant from a random selection
of contestants is chosen for mating.  Typically, each successive generation is
more fit than its predecessor, i.e., the colony's quality improves over time.

To install this package, download the following packages from their URL address:

using Pkg
Pkg.add(url = "https://github.com/tbreloff/ConcreteAbstractions.jl")
Pkg.add(url = "https://github.com/AlanFreed/GeneticAlgorithms.jl")

A genetic algorithm has an interface of

struct GeneticAlgorithm
    c::Colony
end

with constructor

    ga = GeneticAlgorithm(colony)

Method

    run(ga, verbose)

    Method 'run' runs a solver for this genetic algorithm whose population size
    and number of generations to advance through are determined internally.
    A report is written to file for the user to read.  If verbose is 'true,'
    the default, then a page in the report is written for each generation of
    the colony; otherwise, the report only contains information regarding the
    final generation.
"""
module GeneticAlgorithms

using
    Statistics,
    StatsBase

import
    Printf: @sprintf

export
    # macros copied from https://github.com/tbreloff/ConcreteAbstractions.jl
    @base,
    @extend,

    # abstract type
    AbstractSpecies,

    # types, including internal constructors
    Expression,
    Gene,
    Chromosome,
    Genome,
    Creature,
    Colony,
    GeneticAlgorithm,

    # operators
    ==, ≠,

    # external constructors
    procreate,
    alien,
    conceive,

    # methods
    copy,
    deepcopy,
    get,
    getindex,
    set!,
    setindex!,
    toString,

    isDominant,
    isRecessive,
    decode,
    encode!,
    mutate!,
    crossover,
    parameters,
    advanceToNextGeneration!,
    report,
    run,
    runModel,

# constants
    dominant,
    recessive

include("ConcreteAbstractions.jl")

include("Expressions.jl")

include("Genes.jl")

include("Chromosomes.jl")

include("Genomes.jl")

include("Creatures.jl")

include("Species.jl")

include("Colonies.jl")

# The genetic algorithm.

struct GeneticAlgorithm
    c::Colony

    # constructor

    function GeneticAlgorithm(c::Colony)
        new(c)
    end
end # GeneticAlgorithm

# Method

function run(ga::GeneticAlgorithm, verbose::Bool=true)
    println("Each ⋅ represents one generation advanced out of ", ga.c.generationsToConvergence, " generations total.")
    print("    ")

    # Open an IO-stream to write to.
    myDir = string(pwd(), "/files")
    if !isdir(myDir)
        mkdir(myDir)
    end
    myFile = string(myDir, "/ga_report.txt")
    if isfile(myFile)
        myStream = open(myFile; lock=true, read=false, write=true, create=false, truncate=true, append=true)
    else
        myStream = open(myFile; lock=true, read=false, write=true, create=true, truncate=true, append=true)
    end
    seekstart(myStream)

    # The first generation.
    s = string("\n\n", "For generation ", ga.c.generation, " of ", ga.c.generationsToConvergence, ":\n\n")
    s = string(s, report(ga.c))
    write(myStream, s)
    flush(myStream)

    # Run the genetic algorithm.
    if verbose
        for i in 2:ga.c.generationsToConvergence
            print("⋅")
            advanceToNextGeneration!(ga.c)
            s = string("\n\n", "For generation ", ga.c.generation, " of ", ga.c.generationsToConvergence, ":\n\n")
            s = string(s, report(ga.c))
            write(myStream, s)
            flush(myStream)
        end
    else
        for i in 2:ga.c.generationsToConvergence
            print("⋅")
            advanceToNextGeneration!(ga.c)
        end
        s = string("\n\n", "For generation ", ga.c.generation, " of ", ga.c.generationsToConvergence, ":\n\n")
        s = string(s, report(ga.c))
        write(myStream, s)
        flush(myStream)
    end
    close(myStream)
    println()
    println("See ", myFile, " for a report.")
end # run

end # module GeneticAlgorithms
