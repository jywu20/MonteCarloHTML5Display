include("defs.jl")

E_kin_data = Float64[]
mott_data = Float64[]
mag_data = Float64[]

function observe(model::HubbardDQMC, G_up::Matrix{Float64}, G_dn::Matrix{Float64})
    println("observed")
    if double_check
        # Errors between propagated Green functions and Green functions calculated from definition:
        println(relative_err(G_up_τ(model, n_τ), G_up))
        println(relative_err(G_dn_τ(model, n_τ), G_dn))
    end
    
    T_kin = model.T
    n_sites = model.lattice.n_sites

    k_e = - (tr(T_kin * G_up) + tr(T_kin * G_dn)) / n_sites
    push!(E_kin_data, k_e)

    G_up_c = (I - G_up)'
    G_dn_c = (I - G_dn)'
    mott = sum(1 : n_sites) do i
        G_up_c[i, i] * G_dn_c[i, i]
    end
    mott /= n_sites
    push!(mott_data, mott)

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
    push!(mag_data, mag)

    if ! show_binned_only
        @printf "%.10f    %.10f    %.10f \n" k_e mott mag
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
    if show_bin_std 
        println("E_kin      =   $(mean(E_kin_data)) ± $(std(E_kin_data))")
        println("double_occ =   $(mean(mott_data)) ± $(std(mott_data))")
        println("mag        =   $(mean(mag_data)) ± $(std(mag_data))")
    else
        println("E_kin      =   $(mean(E_kin_data)) ")
        println("double_occ =   $(mean(mott_data)) ")
        println("mag        =   $(mean(mag_data)) ")
    end

    println()
end