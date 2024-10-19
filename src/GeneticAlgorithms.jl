#=
Created on Thr 14 Jun 2024
Updated on Fri 18 Oct 2024
Translated from python (code dated 09/08/2017) with enhancements for Julia.
=#

"""
This Julia module implements a genetic algorithm (GA) that was a translation (with refinements and enhancements) of a Python code written by the author for teaching numerical methods to his students at TAMU. This, in turn, was a translation from a Pascal code that, in turn, was a translation of the author's original GA written in the experimental programming language Zonnon. The author used his Zonnon code to illustrate the art of writing object-oriented code to his students at SVSU, for which genetic algorithms are well suited. It turns out that genetic algorithms are also well suited for illustrating the power of Julia's multiple dispatch paradigm. A most notable difference is that the Julia implementation is blazingly fast.

From a biologic interpretation, the genetic algorithm implemented here is a colony of creatures that advances their quality of life (fitness) from one generation to the next. Each instance of type *GeneticAlgorithm* is a colony or collection of creatures whose population is sustained from one generation to the next. Mating between creatures occurs through a process known as tournament play, where the most fit contestant from a random selection of contestants is chosen for mating. Typically, each successive generation is more fit than its predecessor, i.e., the colony's quality improves with each new generation.

To install this package, download the following packages from their URL address:

```julia
using Pkg
Pkg.add(url = "https://github.com/AlanFreed/PhysicalFields.jl")
Pkg.add(url = "https://github.com/AlanFreed/GeneticAlgorithms.jl")
```

The following *GeneticAlgorithm* types can be read-from and written-to a json file.

> *Gene*, *Chromosome*, *Genome*, *Creature*, *Colony*, *ExperimentalData* and *GeneticAlgorithm*.

Users must make their `MyParameters <: AbstractParameters` struct JSON3 compatible in order to be able to read-from and write-to a file. See the examples on how to do this.

Documentation is available for the module *GeneticAlgorithms* (you are reading it), for the exported types (*Gene*, *Chromosome*, *Genome*, *Creature*, *Colony*, *TheData*, *Model*, *AbstractParameters* and *GeneticAlgorithm*), and for the exported constants *dominant* and *recessive*. A template is included in the documentation of *Model*. Documentation for exported methods is included in their type's documentation.

Examples are included in the subdirectory /test.

There is a *README.md* Markdown document that accompanies this software, too.

**Note**: An attempt to *thread* the solver in `_evaluate` lead to bogus results. Maybe someone can figure out how to do this, maybe not.
"""
module GeneticAlgorithms

using
    JSON3,
    PhysicalFields,
    Statistics,
    StatsBase,
    StructTypes,
    Xicor

import
    PhysicalFields as PF

export
    # abstract type
    AbstractParameters,

    # genetic types
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
    parameters,

    # constants
    dominant,
    recessive

include("Genes.jl")

include("Chromosomes.jl")

include("Genomes.jl")

include("Creatures.jl")

include("UserInterface.jl")

include("Colonies.jl")

# The genetic algorithm.

"""
Genetic algorithms are procedures whereby one can parameterize a model against data. Genetic algorithms are probabolistic, in contrast with gradient methods that are deterministic.

To simplify this help, we use the alias
```Julia
import
    PhysicalFields as PF
```

# GeneticAlgorithm

```Julia
struct GeneticAlgorithm
    colony::Colony
end
```
where a `colony` is a population of creatures that evolves from one generation to the next.

## Constructor

```Julia
ga = GeneticAlgorithm(colony)
```

## Method

```Julia
run(ga::GeneticAlgorithm, verbose::Bool = true)
```
Function `run` calls a solver for genetic algorithm `ga` whose population size and number of generations to advance through are determined internally. A report is written to file for the user to read. If `verbose` is `true`, the default, then a page in this report is written for each generation of the colony; otherwise, the report only contains information pertaining to the initial and final generations.

### Persistence

To open or close an IOStream attached to a JSON file, call
```julia
json_stream = PF.openJSONWriter(<my_dir_path>, <my_file_name.json>)
```
> which opens a `json_stream` of type *IOStream* for a file `<my_file_name.json>` located in directory `<my_dir_path>`, both of which are strings, while
```julia
PF.closeJSONStream(json_stream)
```
> flushes the buffer and closes this `json_stream`.

To write or read an instance of type *GeneticAlgorithm* to or from a JSON file, call
```julia
toFile(geneticAlgorithm, json_stream)
```
> which writes a geneticAlgorithm of type *GeneticAlgorithm* to the JSON file attached to a `json_stream` of type *IOStream*, while
```julia
geneticAlgorithm = fromFile(GeneticAlgorithm, json_stream)
```
> reads a geneticAlgorithm of type *GeneticAlgorithm* from the JSON file attached to `json_stream.
"""
struct GeneticAlgorithm
    colony::Colony

    # constructor

    function GeneticAlgorithm(colony::Colony)
        new(colony)::GeneticAlgorithm
    end
end # GeneticAlgorithm

# Methods for storing and retrieving a GeneticAlgorithm to and from a file.

StructTypes.StructType(::Type{GeneticAlgorithm}) = StructTypes.Struct()

function toFile(genetic_algorithm::GeneticAlgorithm, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, genetic_algorithm)
        write(json_stream, '\n')
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    flush(json_stream)
    return nothing
end

function fromFile(::Type{GeneticAlgorithm}, json_stream::IOStream)::GeneticAlgorithm
    if isopen(json_stream)
        genetic_algorithm = JSON3.read(readline(json_stream), GeneticAlgorithm)
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    return genetic_algorithm
end

# Script for execution.

function run!(genetic_algorithm::GeneticAlgorithm; verbose::Bool=true)
    print("Each ⋅ represents one generation evolved out of ")
    print(string(genetic_algorithm.colony.generations))
    println(" generations total.")
    print("  ")

    # Open an IO-stream to write to.
    mydir = string(pwd(), "/files/")
    if !isdir(mydir)
        mkdir(mydir)
    end
    myfile = string(mydir, "ga_report.txt")
    if isfile(myfile)
        mystream = open(myfile; lock=true, read=false, write=true, 
                        create=false, truncate=true, append=true)
    else
        mystream = open(myfile; lock=true, read=false, write=true, 
                        create=true, truncate=true, append=true)
    end
    seekstart(mystream)

    # The first generation.
    s = ""
    g = PF.toString(genetic_algorithm.colony.generation)
    G = string(genetic_algorithm.colony.generations)
    s = string(s, "For generation ", g, " of ", G, ":\n\n")
    s = string(s, report(genetic_algorithm.colony))
    write(mystream, s)
    flush(mystream)

    # Run the genetic algorithm.
    for i in 2:genetic_algorithm.colony.generations-1
        print("⋅")
        advance!(genetic_algorithm.colony)
        if verbose
            g = PF.toString(genetic_algorithm.colony.generation)
            s = string("For generation ", g, " of ", G, ":\n\n")
            s = string(s, report(genetic_algorithm.colony))
            write(mystream, s)
            flush(mystream)
        end
    end
    println("⋅")
    advance!(genetic_algorithm.colony)
    println("The genetic algorithm has finished.\n")
    g = PF.toString(genetic_algorithm.colony.generation)
    s = string("For generation ", g, " of ", G, ":\n\n")
    s = string(s, report(genetic_algorithm.colony))
    write(mystream, s)
    flush(mystream)
    close(mystream)
    println("See file ", myfile, " for a report.")

    return nothing
end # run!

end # module GeneticAlgorithms
