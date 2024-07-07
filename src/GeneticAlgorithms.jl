#=
Created on Thr 14 Jun 2024
Updated on Sun 07 Jul 2024
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
    the colony; otherwise, only a report for the final generation is written.
"""
module GeneticAlgorithms

using
    Statistics,
    StatsBase

import
    Printf: @sprintf

export
    # abstract type
    AbstractSpecies,

    # types, including internal constructors
    Expression,
    Gene,
    Chromosome,
    Genome,
    Creature,
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

    # constants
    dominant,
    recessive

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

# local method

function _writeMyReport(report::String)
    my_dir = string(pwd(), "/files")
    if !isdir(my_dir)
        mkdir(my_dir)
    end
    my_file = string(my_dir, "/ga_report.txt")
    if isfile(my_file)
        my_stream = open(my_file; lock=true, read=false, write=true, create=false, truncate=true, append=true)
    else
        my_stream = open(my_file; lock=true, read=false, write=true, create=true, truncate=true, append=true)
    end
    seekstart(my_stream)
    write(my_stream, report)
    flush(my_stream)
    close(my_stream)
    return nothing
end # _writeMyReport

# Method

function run(ga:GeneticAlgorithm, verbose::Bool=true)
    println("Each ⋅ represents a generation advanced.")
    print("    ")

    # Open an IO-stream to write to.
    my_dir = string(pwd(), "/files")
    if !isdir(my_dir)
        mkdir(my_dir)
    end
    my_file = string(my_dir, "/ga_report.txt")
    if isfile(my_file)
        my_stream = open(my_file; lock=true, read=false, write=true, create=false, truncate=true, append=true)
    else
        my_stream = open(my_file; lock=true, read=false, write=true, create=true, truncate=true, append=true)
    end
    seekstart(my_stream)

    # Run the genetic algorithm.
    if verbose
        s = string("For generation ", ga.c.generation, "of ", ga.c.generationsToConvergence, ":\n\n")
        write(my_stream, s)
        for i in 2:ga.c.generationsToConvergence
            print("⋅")
            advanceToNextGeneration!(ga.c)
            s = string("\n\n", "For generation ", ga.c.generation, "of ", ga.c.generationsToConvergence, ":\n\n")
            s = string(s, report(ga.c))
            write(my_stream, s)
        end
    else
        for i in 2:ga.c.generationsToConvergence-1
            print("⋅")
            advanceToNextGeneration!(ga.c)
        end
        advanceToNextGeneration!(ga.c)
        s = report(ga.c)
        write(my_stream, s)
    end
    println("See ", my_file, " for a report.")
end # run

end # module GeneticAlgorithms
