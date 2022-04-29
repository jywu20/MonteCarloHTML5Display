# Past codes. I need to make it more generalized, for example pertaining to more generic lattices

using LinearAlgebra
using Random
using ProgressMeter
using Statistics

struct HubbardModelDQMC
    U::Float64
    β::Float64
    L::Int64
    Δτ::Float64
    G_up::Array{Float64}
    G_dn::Array{Float64}
    s::Array{Float64}
end

function back_into_range(idx, upper)
    if idx > upper
        return idx % upper
    end
    (idx - upper) % upper + upper
end

function hubbard_dqmc(U::Float64, β::Float64, L::Int64, Δτ::Float64, n_bin::Int64, n_sweep::Int64, n_wrap::Int64, observe)

    # The total number of sites
    n_sites = L^2
    # The total number of imaginary steps
    n_imtimes::Int64 = round(β / Δτ)
    # The decoupling parameter between the auxiliary field and the free fermion field
    α = acosh(exp(Δτ * U / 2))

    # The auxiliary field
    s = rand((-1, 1), (n_imtimes, n_sites))

    inverse_list = zeros(Int64, (L, L))
    site_list = zeros(Int64, (n_sites, 2))
    id = Matrix{Float64}(I, n_sites, n_sites)

    for i in 1:L
        for j in 1:L
            site_list[(i-1)*L+j, :] = [i, j]
            inverse_list[i, j] = (i - 1) * L + j
        end
    end

    neighbor_list = zeros(Int64, (n_sites, 8))
    for i in 1:L
        for j in 1:L
            neighbor_list[inverse_list[i, j], 1] = inverse_list[i, back_into_range(j+1, L)]
            neighbor_list[inverse_list[i, j], 2] = inverse_list[back_into_range(i+1, L), j]
            neighbor_list[inverse_list[i, j], 3] = inverse_list[i, back_into_range(j-1, L)]
            neighbor_list[inverse_list[i, j], 4] = inverse_list[back_into_range(i-1, L), j]
            neighbor_list[inverse_list[i, j], 5] = inverse_list[back_into_range(i+1, L), back_into_range(j+1, L)]
            neighbor_list[inverse_list[i, j], 6] = inverse_list[back_into_range(i+1, L), back_into_range(j-1, L)]
            neighbor_list[inverse_list[i, j], 7] = inverse_list[back_into_range(i-1, L), back_into_range(j-1, L)]
            neighbor_list[inverse_list[i, j], 8] = inverse_list[back_into_range(i-1, L), back_into_range(j+1, L)]
        end
    end

    T = zeros(Float64, (n_sites, n_sites))

    for i in 1:n_sites
        for nn in [1, 2, 3, 4]
            T[i, neighbor_list[i, nn]] = -1.0
        end
    end

    plaquette_sites = zeros(Int64, (2, Int64(n_sites / 4)))

    first_class_count = 1
    second_class_count = 1
    for i in 1:L
        for j in 1:L
            if i % 2 == 1 && j % 2 == 1
                plaquette_sites[1, first_class_count] = inverse_list[i, j]
                first_class_count += 1
            elseif i % 2 == 0 && j % 2 == 0
                plaquette_sites[2, second_class_count] = inverse_list[i, j]
                second_class_count += 1
            end
        end
    end

    function exp_two_hot_sym_mat(val, pos)
        result = copy(id)
        i, j = pos
        result[i, i] = cosh(val)
        result[j, j] = cosh(val)
        result[i, j] = sinh(val)
        result[j, i] = sinh(val)
        result
    end

    exp_kinetic_mat = I

    for m in [1, 2]
        for n in 1:Int64(n_sites/4)
            i1 = plaquette_sites[m, n]
            i2 = neighbor_list[i1, 1]
            i3 = neighbor_list[i1, 5]
            i4 = neighbor_list[i1, 2]
            exp_kinetic_mat = exp_kinetic_mat * exp_two_hot_sym_mat(Δτ, (i1, i2))
            exp_kinetic_mat = exp_kinetic_mat * exp_two_hot_sym_mat(Δτ, (i2, i3))
            exp_kinetic_mat = exp_kinetic_mat * exp_two_hot_sym_mat(Δτ, (i3, i4))
            exp_kinetic_mat = exp_kinetic_mat * exp_two_hot_sym_mat(Δτ, (i4, i1))
        end
    end

    function exp_int_up_mat(s, n)
        diagm(exp.(α * s[n, :]))
    end

    function exp_int_down_mat(s, n)
        diagm(exp.(- α * s[n, :]))
    end

    function b_mat_up_at_tau(s, n)
        exp_int_up_mat(s, n) * exp_kinetic_mat
    end

    function b_mat_down_at_tau(s, n)
        exp_int_down_mat(s, n) * exp_kinetic_mat
    end

    function b_mat_up(s, n2, n1)
        result = id
        for n in (n1+1):n2
            result = b_mat_up_at_tau(s, n) * result
        end
        result
    end

    function b_mat_down(s, n2, n1)
        result = id
        for n in (n1+1):n2
            result = b_mat_down_at_tau(s, n) * result
        end
        result
    end

    function b_mat_up_left(s, n)
        if n == 0
            return (copy(id), copy(id), copy(id))
        end
        U, D, Vt = svd(b_mat_up_at_tau(s, 1))
        D = diagm(D)
        V = Vt'
        for i in 2:n
            U, D, Vt = svd(b_mat_up_at_tau(s, i) * U * D)
            D = diagm(D)
            V = Vt' * V
        end
        (U, D, V)
    end

    function b_mat_up_right(s, n)
        if n == n_imtimes
            return (copy(id), copy(id), copy(id))
        end
        V, D, Ut = svd(b_mat_up_at_tau(s, n_imtimes))
        D = diagm(D)
        U = Ut'
        for i in (n_imtimes-1):-1:(n+1)
            Vt, D, Ut = svd(D * U * b_mat_up_at_tau(s, i))
            U = Ut'
            D = diagm(D)
            V = V * Vt
        end
        (V, D, U)
    end

    function b_mat_down_left(s, n)
        if n == 0
            return (copy(id), copy(id), copy(id))
        end
        U, D, Vt = svd(b_mat_down_at_tau(s, 1))
        D = diagm(D)
        V = Vt'
        for i in 2:n
            U, D, Vt = svd(b_mat_down_at_tau(s, i) * U * D)
            D = diagm(D)
            V = Vt' * V
        end
        (U, D, V)
    end

    function b_mat_down_right(s, n)
        if n == n_imtimes
            return (copy(id), copy(id), copy(id))
        end
        V, D, Ut = svd(b_mat_down_at_tau(s, n_imtimes))
        D = diagm(D)
        U = Ut'
        for i in (n_imtimes-1):-1:(n+1)
            Vt, D, Ut = svd(D * U * b_mat_down_at_tau(s, i))
            U = Ut'
            D = diagm(D)
            V = V * Vt
        end
        (V, D, U)
    end

    function G_up_stable(s, n)
        U_L, D_L, V_L = b_mat_up_left(s, n)
        V_R, D_R, U_R = b_mat_up_right(s, n)
        U, D, Vt = svd(inv(U_R * U_L) + D_L * (V_L * V_R) * D_R)
        V = Vt'
        inv(V * U_R) * diagm(1.0 ./ D) * inv(U_L * U)
    end

    function G_dn_stable(s, n)
        U_L, D_L, V_L = b_mat_down_left(s, n)
        V_R, D_R, U_R = b_mat_down_right(s, n)
        U, D, Vt = svd(inv(U_R * U_L) + D_L * (V_L * V_R) * D_R)
        V = Vt'
        inv(V * U_R) * diagm(1.0 ./ D) * inv(U_L * U)
    end

    function Δ_up(s, n, i)
        exp(-2.0 * α * s[n, i]) - 1.0
    end

    function Δ_dn(s, n, i)
        exp(2.0 * α * s[n, i]) - 1.0
    end

    function accept_ratio_up(s, n, i, G)
        1.0 + Δ_up(s, n, i) * (1.0 - G[i, i])
    end

    function accept_ratio_down(s, n, i, G)
        1.0 + Δ_dn(s, n, i) * (1.0 - G[i, i])
    end

    function G_up_update(s, n, i, G)
        a = G[:, i]
        b = ((I - G)[i, :])'
        G - Δ_up(s, n, i) * a * b / accept_ratio_up(s, n, i, G)
    end

    function G_dn_update(s, n, i, G)
        a = G[:, i]
        b = ((I - G)[i, :])'
        G - Δ_dn(s, n, i) * a * b / accept_ratio_down(s, n, i, G)
    end

    current_G_up = G_up_stable(s, 1)
    current_G_dn = G_dn_stable(s, 1)

    function sweep(tau)
        for i in 1:n_sites

            accept_rate_up = accept_ratio_up(s, tau, i, current_G_up)
            accept_rate_down = accept_ratio_down(s, tau, i, current_G_dn)
            accept_rate = accept_rate_up * accept_rate_down
            
            if rand(Float64) < accept_rate
                current_G_up = G_up_update(s, tau, i, current_G_up)
                current_G_dn = G_dn_update(s, tau, i, current_G_dn)
                s[tau, i] *= -1
            end
        end
    end

    bin_data = []

    for bin_count in 1:n_bin

        sweep_data = []

        for sweep_count in 1:n_sweep
            for tau in 1:(n_imtimes-1)
                sweep(tau)
                
                if tau % n_wrap == 0
                    current_G_up = G_dn_stable(s, tau + 1)
                    current_G_dn = G_dn_stable(s, tau + 1)
                else
                    current_G_up = b_mat_up_at_tau(s, tau + 1) * current_G_up * inv(b_mat_up_at_tau(s, tau + 1))
                    current_G_dn = b_mat_down_at_tau(s, tau + 1) * current_G_dn * inv(b_mat_down_at_tau(s, tau + 1))
                end

            end

            for tau in n_imtimes:-1:2
                sweep(tau)

                if tau % n_wrap == 0
                    current_G_up = G_up_stable(s, tau - 1)
                    current_G_dn = G_dn_stable(s, tau - 1)
                else
                    current_G_up = inv(b_mat_up_at_tau(s, tau)) * current_G_up * b_mat_up_at_tau(s, tau)
                    current_G_dn = inv(b_mat_down_at_tau(s, tau)) * current_G_dn * b_mat_down_at_tau(s, tau)
                end

            end

            push!(sweep_data, observe(current_G_up, current_G_dn, site_list, inverse_list, T))
        end

        push!(bin_data, mean(sweep_data))
    end

    # (
    #     field=s,
    #     current_green_function_up=current_green_function_up,
    #     current_green_function_down=current_green_function_down,
    #     green_function_up_stable=green_function_up_stable,
    #     green_function_down_stable=green_function_down_stable,
    #     green_function_up_update=green_function_up_update,
    #     green_function_down_update=green_function_down_update,
    #     grid_list=grid_list,
    #     inverse_list=inverse_list,
    #     neighbor_list=neighbor_list,
    #     kinetic_hamiltonian = T,
    #     b_mat_up=b_mat_up,
    #     b_mat_down=b_mat_down,
    #     b_mat_up_left=b_mat_up_left,
    #     b_mat_up_right=b_mat_up_right,
    #     b_mat_down_left=b_mat_down_left,
    #     b_mat_down_right=b_mat_down_right,
    #     b_mat_up_at_tau=b_mat_up_at_tau,
    #     b_mat_down_at_tau=b_mat_down_at_tau
    # )
    (grid_list=site_list, inverse_list=inverse_list, neighbor_list=neighbor_list, observables=bin_data)
end

U = 2.0
beta = 4.0
L = 4
dtau = 0.05
n_bin = 10
n_sweep = 400
n_wrap = 10

n_site = L^2
n_imtimes = Int(beta / dtau)

progress = Progress(n_sweep * n_bin)

function observe(G_up, G_down, grid_list, inverse_list, T)
    next!(progress)
    
    E_kin = (tr((I - G_up) * T) + tr((I - G_down) * T)) / n_site

    mott = 0.0
    for i in 1:n_site
        mott += (1.0 - G_up[i, i]) * (1.0 - G_down[i, i]) / n_site
    end
    [E_kin, mott]
end

function relative_err(m1, m2)
    square = x -> x^2
    sum(square.(m1-m2)) / sum(square.(m1))
end

result = hubbard_dqmc(U, beta, L, dtau, n_bin, n_sweep, n_wrap, observe)
observables = reshape(vcat(result.observables...), (2, n_bin))
E_kin_history = observables[1, :]
mott_history = observables[2, :]
println("T $(mean(E_kin_history)) $(std(E_kin_history))")
println("Mott $(mean(mott_history)) $(std(mott_history))")