using LinearAlgebra
using ProgressMeter
using Statistics
using Printf

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

println()
println("Heating up for $heating_up_steps steps.")
sweep!(model, heating_up_steps)
println("Heating up completed.")
println()

#endregion

#region Sweeping

function observable_table_head()
    println("========================================================")
    println("E_kin            double_occ      magnetization")
    println("========================================================")
end

include("observe.jl")

for bin_count in 1 : n_bin
    println("Observables:")
    if ! show_binned_only
        observable_table_head()
    end
    
    sweep!(model, n_sweep; observe = observe)
    println("Bin $bin_count completed.")
    binning()
end

println()
println()
println("Bin average data:")
observable_table_head()
for bin_count in 1 : n_bin
    k_e = E_kin_bin[bin_count]
    mott = mott_bin[bin_count]
    mag = mag_bin[bin_count]
    @printf "%.10f    %.10f    %.10f \n" k_e mott mag
end

println()
println("Binning results:")

println("-------------------------------------------------------------")
println("E_kin      =   $(mean(E_kin_bin)) ± $(std(E_kin_bin))")
println("double_occ =   $(mean(mott_bin)) ± $(std(mott_bin))")
println("mag        =   $(mean(mag_bin)) ± $(std(mag_bin))")

#endregion