#=
Created on Thr 14 Jun 2024
Updated on Wed 14 Aug 2024
Translated from python (code dated 09/08/2017) with enhancements for julia.
=#

"""
This Julia module is a translation of a python module written for the author's
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

To install this package, download the following packages from their URL address:

using Pkg
Pkg.add(url = "https://github.com/AlanFreed/PhysicalFields.jl")
Pkg.add(url = "https://github.com/AlanFreed/GeneticAlgorithms.jl")

JSON files, as implemented by JSON3.jl in Julia, handle the core types:
    Object, Array, String, Number, Bool and Null
where a Number is either a 64-bit integer of a 64-bit floating point number;
hence, all integer fields and all real fields must be described via 64 bits.
A JSON3.Object is an immutable Dict type, while a JSON3.Array is an immutable
Vector type. Command copy(JSON3.Object) will return a mutable Dict object,
while command copy(JSON3.Array) will return a mutable Vector object.

The following GeneticAlgorithm types can be read-from & written-to a json file.
    Gene, Chromosome, Genome, Creature, ExperimentalData, Colony and
    GeneticAlgorithm.
Users must make their MyParameters <: AbstractParameters JSON3 compatible to be
able to read-from and write-to a file.  See the examples on how to do this.

The data structure for genetic algorithm has an interface of

struct GeneticAlgorithm
    c::Colony
end

with constructor

    ga = GeneticAlgorithm(colony)

and procedure

    run(ga, verbose)

    Function 'run' calls a solver for genetic algorithm 'ga' whose population
    size and number of generations to advance through are determined internally.
    A report is written to file for the user to read.  If 'verbose' is true,
    the default, then a page in this report is written for each generation of
    the colony; otherwise, the report only contains information pertaining to
    the final generation.
"""
module GeneticAlgorithms

using
    JSON3,
    PhysicalFields,
    Statistics,
    StatsBase,
    StructTypes

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

    # higher-level methods
    solve,
    advance_to_next_generation!,
    report,
    run,

    # lower-level methods
    copy,
    get,
    getindex,
    set!,
    setindex!,
    toString,
    toBinaryString,
    toFile,
    fromFile,

    isDominant,
    isRecessive,
    mutate!,
    decode,
    encode!,
    crossover,
    phenotypes,

    # constants
    dominant,
    recessive,
    fitness_types

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
        new(colony)::GeneticAlgorithm
    end
end # GeneticAlgorithm

# Methods for storing and retrieving a GeneticAlgorithm to and from a file.

StructTypes.StructType(::Type{GeneticAlgorithm}) = StructTypes.Struct()

"""
Method:\n
    toFile(ga::GeneticAlgorithms.GeneticAlgorithm, json_stream::IOStream)\n
writes a data structure `ga` to the IOStream `json_stream.`\n
For example, consider the code fragment:\n
    json_stream = PhysicalFields.openJSONWriter(<my_dir_path>::String, <my_file_name>::String)\n
    ...\n
    GeneticAlgorithms.toFile(ga::GeneticAlgorithms.GeneticAlgorithm, json_stream::IOStream)\n
    ...\n
    PhysicalFields.closeJSONStream(json_stream::IOStream)\n
where <my_dir_path> is the path to your working directory wherein the file\n
<my_file_name> that is to be written to either exists or will be created,\n
and which must have a .json extension.
"""
function toFile(ga::GeneticAlgorithm, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, ga)
        write(json_stream, '\n')
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    flush(json_stream)
    return nothing
end

"""
Method:\n
    fromFile(::GeneticAlgorithms.GeneticAlgorithm, json_stream::IOStream)\n
reads an instance of type GeneticAlgorithm from the IOStream `json_stream.`\n
For example, consider the code fragment:\n
    json_stream = PhysicalFields.openJSONReader(<my_dir_path>::String, <my_file_name>::String)\n
    ...\n
    ga = GeneticAlgorithms.fromFile(::GeneticAlgorithms.GeneticAlgorithm, json_stream::IOStream)\n
    ...\n
    PhysicalFields.closeJSONStream(json_stream::IOStream)\n
which returns a `ga,` an object of type GeneticAlgorithms.GeneticAlgorithm.\n
Here <my_dir_path> is the path to your working directory wherein the file\n
to be read from, i.e., <my_file_name>, must exist, and which is to have a\n
.json extension.
"""
function fromFile(::Type{GeneticAlgorithm}, json_stream::IOStream)::GeneticAlgorithm
    if isopen(json_stream)
        ga = JSON3.read(readline(json_stream), GeneticAlgorithm)
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    return ga
end

# Script for execution.

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
