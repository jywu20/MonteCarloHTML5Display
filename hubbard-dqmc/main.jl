# The main Monte Carlo loop
using Statistics
using Printf

include("lib.jl")

#region Model parameters

t = 1.0
U = 4.0
β = 4.0

#endregion

#region Problem size

L = 10 
lattice = SquareLattice2D(L)
n_sites = lattice.n_sites

n_imtimes = 100
n_wrap = 10

#endregion

#region Sweeping control

n_heat = 1000
n_sweep = 4000
n_bin = 10

#endregion

E_kin_data = zeros(n_bin)
mott_data = zeros(n_bin)

model = HubbardDQMC(SquareLattice2D, lattice, t, U, β, n_imtimes, n_wrap)
T = model.T

##

println("Square lattice 2D Hubbard simple DQMC")
println()
println()

println("t   =    $t")
println("U   =    $U")
println("β   =    $β")
println()
println("n_sites   =   $n_sites")
println("n_imtimes =   $n_imtimes")
println("n_wrap    =   $n_wrap")
println()
println("n_heat    =   $n_heat")
println("n_sweep   =   $n_sweep")
println("n_bi n    =   $n_bin")

println()

sweep!(model, n_heat)

println("Heating up finished.")
println()
println()

println("E_kin         double_occ")
println("====================================")

for bin_count in 1 : n_bin
    E_kin_history = zeros(n_sweep)
    mott_history = zeros(n_sweep)

    sweep_count = 0
    sweep!(model, n_sweep) do G_up, G_dn, model
        sweep_count += 1

        E_kin = (tr((I - G_up) * T) + tr((I - G_dn) * T)) / n_sites

        mott = 0.0
        for i in 1:n_sites
            mott += (1.0 - G_up[i, i]) * (1.0 - G_dn[i, i]) / n_sites
        end

        E_kin_history[sweep_count] = E_kin
        mott_history[sweep_count] = mott
    end

    E_kin_data[bin_count] = mean(E_kin_history)
    mott_data[bin_count] = mean(mott_history)

    @printf "%.10f    %.10f \n" E_kin_data[bin_count] mott_data[bin_count] 
end

println()
println("====================================")
println()

println("E_kin      =   $(mean(E_kin_data)) ± $(std(E_kin_data))")
println("double_occ =   $(mean(mott_data)) ± $(std(mott_data))")