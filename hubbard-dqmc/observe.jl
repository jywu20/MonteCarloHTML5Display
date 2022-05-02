include("defs.jl")

observables_names = ["Kinetic energy"]
observables_history = zeros(length(observables_names), n_bin)

kinetic_energy_history = Float64[]
mott_history = Float64[]
mag_history = Float64[]

function observe()
    if double_check
        # Errors between propagated Green functions and Green functions calculated from definition:
        open(error_path, "a") do error_file
            println(error_file, relative_err(G_up_τ(τ_now), G_up))
            println(error_file, relative_err(G_dn_τ(τ_now), G_dn))
        end
    end
    
    open(observables_path, "a") do observable_file
        k_e = - (tr(T_kin * G_up) + tr(T_kin * G_dn)) / n_sites
        push!(kinetic_energy_history, k_e)

        G_up_c = (I - G_up)'
        G_dn_c = (I - G_dn)'
        mott = sum(1 : n_sites) do i
            G_up_c[i, i] * G_dn_c[i, i]
        end
        mott /= n_sites
        push!(mott_history, mott)

        mag = sum(Iterators.product(1 : n_sites, 1 : n_sites)) do point_pair
            i, j = point_pair
            r_i = lattice.site_list[i, :]
            r_j = lattice.site_list[j, :]
            fourier_prefactor = (-1)^sum(r_i - r_j)
            sᶻ_isᶻ_j = (G_up_c[i, i] * G_up_c[j, j] + G_up_c[i, j] * G_up[i, j] 
                      + G_dn_c[i, i] * G_dn_c[j, j] + G_dn_c[i, j] * G_dn[i, j]
                      - G_up_c[i, i] * G_dn_c[j, j] 
                      - G_dn_c[i, i] * G_up_c[j, j])
            fourier_prefactor * sᶻ_isᶻ_j
        end
        # The factor 4 comes from spin 1/2
        mag /= 4 * n_sites^2
        push!(mag_history, mag)

        if ! binned_only
            println(observable_file, "$k_e    $mott    $mag")
        end
    end
end

function mean_no_nan(array::Vector{Float64})::Float64
    if any(isnan, array)
        dqmc_log("Warning: NaN occurs in binning.")
        return mean(filter(x -> ! isnan(x), array))
    end
    return mean(array)
end

function filter_out_nan(array::Vector{Float64})::Float64
    filter(x -> ! isnan(x), array)
end

function autocorrelation(seq::Vector{Float64}, span::Int64)::Float64
    m = mean(seq)
    sum((seq[span + 1 : end] .- m) .* (seq[1 : end - span] .- m)) / (length(seq) - span)
end

function binning()
    open(observables_path, "a") do observable_file
        if show_bin_std 
            println(observable_file, 
            "Binned Kinetic energy: $(mean_no_nan(kinetic_energy_history)) ± $(std(kinetic_energy_history))")
            println(observable_file, 
            "Binned double occupation: $(mean_no_nan(mott_history)) ± $(std(mott_history))")
            println(observable_file, 
            "Magnetization: $(mean_no_nan(mag_history)) ± $(std(mag_history))")
        else
            println(observable_file, "Binned Kinetic energy: $(mean_no_nan(kinetic_energy_history))")
        end
    end
end