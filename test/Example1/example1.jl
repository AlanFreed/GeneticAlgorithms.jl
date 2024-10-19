module example1

# To use this code, you will need to download the genetic algorithm at:
# using Pkg
# Pkg.add(url = "https://github.com/AlanFreed/PhysicalFields.jl")
# Pkg.add(url = "https://github.com/AlanFreed/GeneticAlgorithms.jl")

using
    CairoMakie,
    GeneticAlgorithms,
    PhysicalFields

import
    GeneticAlgorithms as GA,
    PhysicalFields    as PF

export
    run
#=
--------------------------------------------------------------------------------
=#

function myData()::GA.TheData
    n_exp = 1
    pairs = 1
    n_pts = 20
    
    experiments     = n_exp
    conjugate_pairs = Vector{Int}(undef, n_exp)
    data_points     = Vector{Int}(undef, n_exp)
    for exp in 1:n_exp
        conjugate_pairs[exp] = pairs
        data_points[exp]     = n_pts
    end
    
    time     = Vector{PF.ArrayOfPhysicalScalars}(undef, 0) # a static model
    control  = Vector{Vector{PF.ArrayOfPhysicalScalars}}(undef, n_exp)
    response = Vector{Vector{PF.ArrayOfPhysicalScalars}}(undef, n_exp)
    for exp in 1:n_exp
        control[exp]  = Vector{PF.ArrayOfPhysicalScalars}(undef, pairs)
        response[exp] = Vector{PF.ArrayOfPhysicalScalars}(undef, pairs)
        for pair in 1:pairs
            control[exp][pair]  = PF.ArrayOfPhysicalScalars(n_pts, PF.STRAIN)
            response[exp][pair] = PF.ArrayOfPhysicalScalars(n_pts, PF.STRESS)
        end
    end
    control[1][1][1]  = PF.PhysicalScalar(0.1, PF.STRAIN)
    control[1][1][2]  = PF.PhysicalScalar(0.2, PF.STRAIN)
    control[1][1][3]  = PF.PhysicalScalar(0.3, PF.STRAIN)
    control[1][1][4]  = PF.PhysicalScalar(0.4, PF.STRAIN)
    control[1][1][5]  = PF.PhysicalScalar(0.5, PF.STRAIN)
    control[1][1][6]  = PF.PhysicalScalar(0.6, PF.STRAIN)
    control[1][1][7]  = PF.PhysicalScalar(0.7, PF.STRAIN)
    control[1][1][8]  = PF.PhysicalScalar(0.8, PF.STRAIN)
    control[1][1][9]  = PF.PhysicalScalar(0.9, PF.STRAIN)
    control[1][1][10] = PF.PhysicalScalar(1.0, PF.STRAIN)
    control[1][1][11] = PF.PhysicalScalar(1.1, PF.STRAIN)
    control[1][1][12] = PF.PhysicalScalar(1.2, PF.STRAIN)
    control[1][1][13] = PF.PhysicalScalar(1.3, PF.STRAIN)
    control[1][1][14] = PF.PhysicalScalar(1.4, PF.STRAIN)
    control[1][1][15] = PF.PhysicalScalar(1.5, PF.STRAIN)
    control[1][1][16] = PF.PhysicalScalar(1.6, PF.STRAIN)
    control[1][1][17] = PF.PhysicalScalar(1.7, PF.STRAIN)
    control[1][1][18] = PF.PhysicalScalar(1.8, PF.STRAIN)
    control[1][1][19] = PF.PhysicalScalar(1.9, PF.STRAIN)
    control[1][1][20] = PF.PhysicalScalar(2.0, PF.STRAIN)
    
    response[1][1][1]  = PF.PhysicalScalar(0.050, PF.STRESS)
    response[1][1][2]  = PF.PhysicalScalar(0.111, PF.STRESS)
    response[1][1][3]  = PF.PhysicalScalar(0.193, PF.STRESS)
    response[1][1][4]  = PF.PhysicalScalar(0.290, PF.STRESS)
    response[1][1][5]  = PF.PhysicalScalar(0.349, PF.STRESS)
    response[1][1][6]  = PF.PhysicalScalar(0.450, PF.STRESS)
    response[1][1][7]  = PF.PhysicalScalar(0.559, PF.STRESS)
    response[1][1][8]  = PF.PhysicalScalar(0.622, PF.STRESS)
    response[1][1][9]  = PF.PhysicalScalar(0.744, PF.STRESS)
    response[1][1][10] = PF.PhysicalScalar(0.835, PF.STRESS)
    response[1][1][11] = PF.PhysicalScalar(1.032, PF.STRESS)
    response[1][1][12] = PF.PhysicalScalar(1.144, PF.STRESS)
    response[1][1][13] = PF.PhysicalScalar(1.266, PF.STRESS)
    response[1][1][14] = PF.PhysicalScalar(1.396, PF.STRESS)
    response[1][1][15] = PF.PhysicalScalar(1.409, PF.STRESS)
    response[1][1][16] = PF.PhysicalScalar(1.494, PF.STRESS)
    response[1][1][17] = PF.PhysicalScalar(1.625, PF.STRESS)
    response[1][1][18] = PF.PhysicalScalar(1.675, PF.STRESS)
    response[1][1][19] = PF.PhysicalScalar(1.700, PF.STRESS)
    response[1][1][20] = PF.PhysicalScalar(1.710, PF.STRESS)
    
    return GA.TheData(experiments, conjugate_pairs, data_points, 
                      time, control, response)
end # myData

mutable struct MyParameters <: GA.AbstractParameters
    a::PF.PhysicalScalar
    b::PF.PhysicalScalar
    E::PF.PhysicalScalar
    
    function MyParameters()
        a = PF.PhysicalScalar(PF.DIMENSIONLESS)
        b = PF.PhysicalScalar(PF.DIMENSIONLESS)
        E = PF.PhysicalScalar(PF.MODULUS)
        new(a, b, E)::MyParameters
    end
end # MyParameters

myparameters = MyParameters()
mydata  = myData()
mymodel = GA.Model{MyParameters}(myparameters::MyParameters,
                                 mydata::GA.TheData)

function _mysolve(mymodel::GA.Model;
                  exp::Int, pair::Int, datum::Int)::PF.PhysicalScalar
    # Example of a static model described by a simple function, i.e., σ = f(ε).
    ε = mymodel.data.control[exp][pair][datum]
    θ = mymodel.θ
    σ = θ.E * sinh(θ.a*ε) / (θ.a * cosh(θ.a*ε) - θ.b * sinh(θ.a*ε))
    return σ
end # _mysolve

function GA.solve!(::MyParameters, mymodel::GA.Model)
    # Type MyParameters allows Julia to select solve! via multiple dispatch.
    for exp in 1:mymodel.data.experiments
        for pair in 1:mymodel.data.conjugate_pairs[exp]
            for datum in 1:mymodel.data.data_points[exp]
                y = _mysolve(mymodel; exp, pair, datum)
                mymodel.data.response_mod[exp][pair][datum] = y
            end
        end
    end
    return nothing
end # GA.solve!

function plot(colony::Colony)
    # Create a figure illustrating a best fit of the model.
    εₑ = Vector{Float64}
    σₑ = Vector{Float64}
    σₘ = Vector{Float64}
    for exp in 1:colony.mydata.experiments
        εₑ = zeros(Float64, colony.mydata.data_points[exp])
        σₑ = zeros(Float64, colony.mydata.data_points[exp])
        σₘ = zeros(Float64, colony.mydata.data_points[exp])
        for pair in 1:colony.mydata.conjugate_pairs[exp]
            for datum in 1:colony.mydata.data_points[exp]
                εₑ[datum] = PF.get(colony.mydata.control[exp][pair][datum])
                σₑ[datum] = PF.get(colony.mydata.response_exp[exp][pair][datum])
                σₘ[datum] = PF.get(colony.mydata.response_mod[exp][pair][datum])
            end
        end
    end
    
    σᵢ    = zeros(Float64, 2)
    σᵢ[2] = PF.get(colony.mydata.response_exp[1][1][colony.mydata.data_points[1]])
    
    fig = Figure(; size = (1000, 500))
    ax1 = Axis(fig[1, 1];
        xlabel = "strain ε",
        ylabel = "stress σ (Pa)",
        title = "Genetic Algorithm",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    scatter!(ax1, εₑ, σₑ;
        marker = :circle,
        markersize = 10,
        color = :black,
        label = "data")
    lines!(ax1, εₑ, σₘ;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "model")
    axislegend("Legend",
        position = :lt)
    ax2 = Axis(fig[1, 2];
        xlabel = "experimental stress σₑ (Pa)",
        ylabel = "model for stress σₘ (Pa)",
        title = "Fit to Data",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    scatter!(ax2, σₑ, σₘ;
        marker = :circle,
        markersize = 10,
        color = :black,
        label = "σₑ vs. σₘ")
    lines!(ax2, σᵢ, σᵢ;
        linewidth = 1,
        linestyle = :solid,
        color = :black,
        label = "perfect fit")
    axislegend("Legend",
        position = :lt)
    save(string(pwd(), "/files/GA_example1.png"), fig)
end # plot

function run()
    # Specify it the report is to be verbose or not.
    verbose = false

    myparameters           = MyParameters()
    mydata                 = myData()
    probability_mutation   = 0.01 
    probability_crossover  = 0.95
    probability_immigrant  = 0.005
    parameters_alien       = Vector{PF.PhysicalScalar}(undef, 0)    # no alien
    parameters_min         = Vector{PF.PhysicalScalar}(undef, 3)
    parameters_max         = Vector{PF.PhysicalScalar}(undef, 3)
    parameters_constrained = Vector{Tuple{Int,Int}}(undef, 1) 
    significant_figures    = 6

    # Populate parameter bounds and constraints.
    parameters_min[1] = PF.PhysicalScalar(0.1, PF.DIMENSIONLESS)
    parameters_min[2] = PF.PhysicalScalar(0.1, PF.DIMENSIONLESS)
    parameters_min[3] = PF.PhysicalScalar(0.1, PF.MODULUS)

    parameters_max[1] = PF.PhysicalScalar(1.5, PF.DIMENSIONLESS)
    parameters_max[2] = PF.PhysicalScalar(1.5, PF.DIMENSIONLESS)
    parameters_max[3] = PF.PhysicalScalar(1.5, PF.MODULUS)

    parameters_constrained[1] = (2, 1)  # b < a
    
    # Create the colony and its genetic algorithm.
    colony = GA.Colony(myparameters, mydata, probability_mutation,
                       probability_crossover, probability_immigrant,
                       parameters_alien, parameters_min, parameters_max,
                       parameters_constrained, significant_figures)

    genetic_algorithm = GA.GeneticAlgorithm(colony)
    GA.run!(genetic_algorithm; verbose)
    plot(genetic_algorithm.colony)
end # run

end # example1
