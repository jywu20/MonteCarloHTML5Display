struct HubbardDQMC{L <: AbstractLattice}
    U::Float64
    β::Float64
    Δτ::Float64
    n_imtimes::Int64
    n_wrap::Int64
    α::Float64
    # G[i, j, τ]
    G_up::Array{Float64, 3}
    G_dn::Array{Float64, 3}
    lattice::L
    # s[τ, x]
    s::Matrix{Int64}
    # T[i, j]
    T::Matrix{Float64}
    id::Matrix{Float64}
    exp_kinetic_mat::Matrix{Float64}
end

function HubbardDQMC(::Type{L}, lattice::L, U::Float64, β::Float64, n_imtimes::Int64, n_wrap::Int64; 
    s = nothing) where {L <: AbstractLattice}
    Δτ = β / n_imtimes
    α = acosh(exp(Δτ * U / 2))

    # The Green functions are initialized as zeros
    n_sites = lattice.n_sites
    neighbor_list = lattice.neighbor_list
    inverse_list = lattice.inverse_list
    site_list = lattice.site_list

    G_up = zeros(n_sites, n_sites, n_imtimes)
    G_dn = zeros(n_sites, n_sites, n_imtimes)
    
    # If no initial configuration is passed into this constructor, we use a random one
    if s === nothing
        s = rand((-1, 1), (n_imtimes, n_sites))
    end

    # Building the kinetic energy matrix
    T = zeros(Float64, (n_sites, n_sites))

    for i in 1:n_sites
        # Nearest hopping 
        for nn in lattice.neighbor_list_indices[1]
            T[i, neighbor_list[i, nn]] = -1.0
        end
    end

    id = Matrix{Float64}(I, n_sites, n_sites)

    plaquette_sites = zeros(Int64, (2, Int64(n_sites / 4)))

    first_class_count = 1
    second_class_count = 1
    
    for site_idx in 1 : n_sites
        i, j = site_list[site_idx, :]
        if i % 2 == 1 && j % 2 == 1
            plaquette_sites[1, first_class_count] = inverse_list[i, j]
            first_class_count += 1
        elseif i % 2 == 0 && j % 2 == 0
            plaquette_sites[2, second_class_count] = inverse_list[i, j]
            second_class_count += 1
        end
    end

    exp_kinetic_mat = I

    for m in [1, 2]
        for n in 1:Int64(n_sites/4)
            i1 = plaquette_sites[m, n]
            i2 = neighbor_list[i1, 1]
            i3 = neighbor_list[i1, 5]
            i4 = neighbor_list[i1, 2]
            exp_kinetic_mat = exp_kinetic_mat * exp_two_hot_sym_mat(Δτ, (i1, i2), id)
            exp_kinetic_mat = exp_kinetic_mat * exp_two_hot_sym_mat(Δτ, (i2, i3), id)
            exp_kinetic_mat = exp_kinetic_mat * exp_two_hot_sym_mat(Δτ, (i3, i4), id)
            exp_kinetic_mat = exp_kinetic_mat * exp_two_hot_sym_mat(Δτ, (i4, i1), id)
        end
    end

    HubbardDQMC(U, β, Δτ, n_imtimes, n_wrap, α, G_up, G_dn, lattice, s, T, id, exp_kinetic_mat)
end