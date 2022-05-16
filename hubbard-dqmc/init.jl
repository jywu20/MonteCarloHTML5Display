function HubbardDQMC(::Type{L}, lattice::L, t::Float64, U::Float64, β::Float64, n_τ::Int64, n_wrap::Int64; 
    s = nothing) where {L <: AbstractLattice}
    Δτ = β / n_τ
    α = acosh(exp(Δτ * U / 2))

    n_sites = lattice.n_sites
    neighbor_list = lattice.neighbor_list
    inverse_list = lattice.inverse_list
    site_list = lattice.site_list

    # If no initial configuration is passed into this constructor, we use a random one
    if s === nothing
        s = rand((-1, 1), (n_τ, n_sites))
    end

    # Building the kinetic energy matrix
    T = zeros(Float64, (n_sites, n_sites))

    for i in 1:n_sites
        # Nearest hopping 
        for nn in lattice.neighbor_list_indices[1]
            T[i, neighbor_list[i, nn]] = - t
        end
    end

    id = Matrix{Float64}(I, n_sites, n_sites)

    #plaquette_sites = zeros(Int64, (2, Int64(n_sites / 4)))
    #
    #first_class_count = 1
    #second_class_count = 1
    #
    #for site_idx in 1 : n_sites
    #    i, j = site_list[site_idx, :]
    #    if i % 2 == 1 && j % 2 == 1
    #        plaquette_sites[1, first_class_count] = inverse_list[i, j]
    #        first_class_count += 1
    #    elseif i % 2 == 0 && j % 2 == 0
    #        plaquette_sites[2, second_class_count] = inverse_list[i, j]
    #        second_class_count += 1
    #    end
    #end

    # exp_kinetic_mat = I
    #for m in [1, 2]
    #    for n in 1:Int64(n_sites/4)
    #        i1 = plaquette_sites[m, n]
    #        i2 = neighbor_list[i1, 1]
    #        i3 = neighbor_list[i1, 5]
    #        i4 = neighbor_list[i1, 2]
    #        exp_kinetic_mat = exp_kinetic_mat * exp_two_hot_sym_mat(Δτ, (i1, i2), id)
    #        exp_kinetic_mat = exp_kinetic_mat * exp_two_hot_sym_mat(Δτ, (i2, i3), id)
    #        exp_kinetic_mat = exp_kinetic_mat * exp_two_hot_sym_mat(Δτ, (i3, i4), id)
    #        exp_kinetic_mat = exp_kinetic_mat * exp_two_hot_sym_mat(Δτ, (i4, i1), id)
    #    end
    #end

    # Stroed B-matrices, to avoid unncessisary calculations.
    B_up_storage = zeros(n_sites, n_sites, n_τ)
    B_dn_storage = zeros(n_sites, n_sites, n_τ)

    model = HubbardDQMC(t, U, β, 
    Δτ, n_τ, 
    n_wrap, 
    α, 
    lattice, s, 
    T, exp(Δτ * T), exp(- Δτ * T), 
    id, 
    B_up_storage, B_dn_storage)

    for τ in 1 : n_τ
        B_up_storage[:, :, τ] = B_up(model, τ)
        B_dn_storage[:, :, τ] = B_dn(model, τ)
    end

    model
end
