# The main Monte Carlo loop
using Statistics
using Printf

include("lib.jl")

#region Model parameters

t = 1.0
U = 8.0
β = 4.0

#endregion

#region Problem size

L = 4 
lattice = SquareLattice2D(L)
n_sites = lattice.n_sites
site_list = lattice.site_list

Δτ = 0.05
n_imtimes = Int(β / Δτ) 
n_wrap = 10

#endregion

#region Sweeping control

n_heat = Int(4 * n_sites * β) 
n_sweep = 500
n_bin = 10

#endregion

E_kin_data = zeros(n_bin)
mott_data = zeros(n_bin)
mag_data = zeros(n_bin)

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
println("n_bin     =   $n_bin")

println()

flush(stdout)

sweep!(model, n_heat)

println("Heating up finished.")
println()
println()

println("====================================================")
println("E_kin            double_occ    magnetization")
println("====================================================")

flush(stdout)

for bin_count in 1 : n_bin
    E_kin_history = zeros(n_sweep)
    mott_history = zeros(n_sweep)
    mag_history = zeros(n_sweep)

    sweep_count = 0
    sweep!(model, n_sweep) do G_up, G_dn, model
        sweep_count += 1
        # ⟨c†_j c_i⟩ = δ_{ij} - ⟨c_i c†_j⟩, and G_up = [⟨c_i c†_j⟩]_{ij}
        G_up_c = (I - G_up)'
        G_dn_c = (I - G_dn)'

        E_kin = (tr(G_up_c * T) + tr(G_dn_c * T)) / n_sites

        mott = sum(1 : n_sites) do i
            G_up_c[i, i] * G_dn_c[i, i]
        end
        mott /= n_sites

        Q = [π, π]
        mag = sum(Iterators.product(1 : n_sites, 1 : n_sites)) do point_pair
            i, j = point_pair
            r_i = site_list[i, :]
            r_j = site_list[j, :]
            fourier_prefactor = exp(- im * Q' * (r_i - r_j))
            sᶻ_isᶻ_j = (G_up_c[i, i] * G_up_c[j, j] + G_up_c[i, j] * G_up[i, j] 
                      + G_dn_c[i, i] * G_dn_c[j, j] + G_dn_c[i, j] * G_dn[i, j]
                      - G_up_c[i, i] * G_dn_c[j, j] 
                      - G_dn_c[i, i] * G_up_c[j, j])
            fourier_prefactor * sᶻ_isᶻ_j
        end
        # The factor 4 comes from spin 1/2
        mag /= 4 * n_sites^2

        E_kin_history[sweep_count] = E_kin
        mott_history[sweep_count] = mott
        # A small imaginary part may occur
        mag_history[sweep_count] = abs(mag)
    end

    E_kin_data[bin_count] = mean(E_kin_history)
    mott_data[bin_count] = mean(mott_history)
    mag_data[bin_count] = mean(mag_history)

    @printf "%.10f    %.10f    %.10f \n" E_kin_data[bin_count] mott_data[bin_count] mag_data[bin_count]
    flush(stdout)
end

println()
println("====================================================")
println()

println("E_kin      =   $(mean(E_kin_data)) ± $(std(E_kin_data))")
println("double_occ =   $(mean(mott_data)) ± $(std(mott_data))")
println("mag        =   $(mean(mag_data)) ± $(std(mag_data))")