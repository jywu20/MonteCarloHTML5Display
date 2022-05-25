using LinearAlgebra
using ProgressMeter
using Statistics
using Printf
using Plots

include("lib.jl")

include("config.jl")

#region Initial logging

println("Welcome to a simple DQMC simulation of Hubbard model.")
println()

lattice = SquareLattice2D(n_side)
model = HubbardDQMC(SquareLattice2D, lattice, t, U, β, n_τ, n_wrap)

println("Initialization completed.")

# Print parameters used in the simulation.
println("We are simulating with")
println("  t      =   $(string(t))")
println("  U      =   $(string(U))")
println("  beta   =   $(string(β))")
println("  dtau   =   $(string(model.Δτ))")
println("  alpha  =   $(string(model.α))")
println("on a $n_side * $n_side 2d lattice and with $n_τ time steps.")

if double_check
    println()
    println("Double checking is on. During the calculation, the program double-checks the error between the optimized version of B-matrices, Green functions, etc. and their by-definition values.")
    println("This can be extremly slow. Turn off the marker when scaling up.")
    println()
end

if double_check
    println("numerical error found in the double-check process")
end

#endregion

#region Heating up

heating_up_k_data = zeros(heating_up_steps)

let 
    progress = Progress(heating_up_steps)
    heating_count = 0

    println()
    println("Heating up for $heating_up_steps steps.")
    sweep!(model, heating_up_steps) do model, G_up, G_dn
        heating_count += 1

        T_kin = model.T
        n_sites = model.lattice.n_sites
        k_e = - (tr(T_kin * G_up) + tr(T_kin * G_dn)) / n_sites
        heating_up_k_data[heating_count] = k_e
        println(k_e)

        #next!(progress)
    end
    println("Heating up completed.")
    println()
end

#endregion

#regoin Show configuration 

working_path = "D:/Projects/Modern Physics Experiments/HTML5/hubbard-report/configuration-u-$(Int(U))/"

##
n_sweep = 100
let 
    sweep_count = 0
    sweep!(model, n_sweep) do model, G_up, G_dn
        sweep_count += 1
        heatmap(reshape(model.s[30, :], (n_side, n_side)), 
        c = cgrad([:orange, :blue]), legend = :none, aspect_ratio=1, 
        xlims = (0.49, n_side + 0.5), ylims = (0.5, n_side + 0.5), dpi = 500)
        savefig(working_path * "fixed-tau-markov-evolution-$sweep_count.png")
    end
end

#endregion

#regoin Show configuration, but this time update τ, no Markov update

working_path = "D:/Projects/Modern Physics Experiments/HTML5/hubbard-report/configuration-u-$(Int(U))/"

##

for τ in 1 : n_τ
    heatmap(reshape(model.s[τ, :], (n_side, n_side)), 
    c = cgrad([:orange, :blue]), legend = :none, aspect_ratio=1, 
    xlims = (0.49, n_side + 0.5), ylims = (0.5, n_side + 0.5), dpi = 500)
    savefig(working_path * "no-markov-evolution-tau-$τ.png")
end

#endregion