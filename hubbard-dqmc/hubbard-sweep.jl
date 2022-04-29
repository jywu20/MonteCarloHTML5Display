function exp_int_up_mat(model::HubbardDQMC, n::Int64)
    s = model.s
    α = model.α
    diagm(exp.(α * s[n, :]))
end

function exp_int_down_mat(model::HubbardDQMC, n::Int64)
    s = model.s
    α = model.α
    diagm(exp.(- α * s[n, :]))
end

function B_up_τ(model::HubbardDQMC, n::Int64)
    exp_kinetic_mat = model.exp_kinetic_mat
    exp_int_up_mat(model, n) * exp_kinetic_mat
end

function B_dn_τ(model::HubbardDQMC, n::Int64)
    exp_kinetic_mat = model.exp_kinetic_mat
    exp_int_down_mat(model, n) * exp_kinetic_mat
end

function B_up(model::HubbardDQMC, n2::Int64, n1::Int64)
    result = model.id
    for n in (n1+1):n2
        result = B_up_τ(model, n) * result
    end
    result
end

function B_dn(model::HubbardDQMC, n2::Int64, n1::Int64)
    result = model.id
    for n in (n1+1):n2
        result = B_dn_τ(model, n) * result
    end
    result
end

function B_up_left(model::HubbardDQMC, n::Int64)
    id = model.id
    if n == 0
        return (copy(id), copy(id), copy(id))
    end
    U, D, Vt = svd(B_up_τ(model, 1))
    D = diagm(D)
    V = Vt'
    for i in 2:n
        U, D, Vt = svd(B_up_τ(model, i) * U * D)
        D = diagm(D)
        V = Vt' * V
    end
    (U, D, V)
end

function B_up_right(model::HubbardDQMC, n::Int64)
    id = model.id
    if n == n_imtimes
        return (copy(id), copy(id), copy(id))
    end
    V, D, Ut = svd(B_up_τ(model, n_imtimes))
    D = diagm(D)
    U = Ut'
    for i in (n_imtimes-1):-1:(n+1)
        Vt, D, Ut = svd(D * U * B_up_τ(model, i))
        U = Ut'
        D = diagm(D)
        V = V * Vt
    end
    (V, D, U)
end

function B_dn_left(model::HubbardDQMC, n::Int64)
    if n == 0
        return (copy(id), copy(id), copy(id))
    end
    U, D, Vt = svd(B_dn_τ(model, 1))
    D = diagm(D)
    V = Vt'
    for i in 2:n
        U, D, Vt = svd(B_dn_τ(model, i) * U * D)
        D = diagm(D)
        V = Vt' * V
    end
    (U, D, V)
end

function B_dn_right(model::HubbardDQMC, n::Int64)
    id = model.id 
    if n == n_imtimes
        return (copy(id), copy(id), copy(id))
    end
    V, D, Ut = svd(B_dn_τ(model, n_imtimes))
    D = diagm(D)
    U = Ut'
    for i in (n_imtimes-1):-1:(n+1)
        Vt, D, Ut = svd(D * U * B_dn_τ(model, i))
        U = Ut'
        D = diagm(D)
        V = V * Vt
    end
    (V, D, U)
end

function G_up_stable(model::HubbardDQMC, n::Int64)
    U_L, D_L, V_L = B_up_left(model, n)
    V_R, D_R, U_R = B_up_right(model, n)
    U, D, Vt = svd(inv(U_R * U_L) + D_L * (V_L * V_R) * D_R)
    V = Vt'
    inv(V * U_R) * diagm(1.0 ./ D) * inv(U_L * U)
end

function G_dn_stable(model::HubbardDQMC, n::Int64)
    U_L, D_L, V_L = B_dn_left(model, n)
    V_R, D_R, U_R = B_dn_right(model, n)
    U, D, Vt = svd(inv(U_R * U_L) + D_L * (V_L * V_R) * D_R)
    V = Vt'
    inv(V * U_R) * diagm(1.0 ./ D) * inv(U_L * U)
end

function Δ_up(model::HubbardDQMC, n::Int64, i::Int64)
    α = model.α
    s = model.s
    exp(-2.0 * α * s[n, i]) - 1.0
end

function Δ_dn(model::HubbardDQMC, n::Int64, i::Int64)
    α = model.α
    s = model.s
    exp(2.0 * α * s[n, i]) - 1.0
end

function accept_ratio_up(model::HubbardDQMC, n::Int64, i::Int64, G::Matrix{Float64})
    1.0 + Δ_up(model, n, i) * (1.0 - G[i, i])
end

function accept_ratio_down(model::HubbardDQMC, n::Int64, i::Int64, G::Matrix{Float64})
    1.0 + Δ_dn(model, n, i) * (1.0 - G[i, i])
end

function G_up_update(model::HubbardDQMC, n::Int64, i::Int64, G::Matrix{Float64})
    a = G[:, i]
    b = ((I - G)[i, :])'
    G - Δ_up(model, n, i) * a * b / accept_ratio_up(model, n, i, G)
end

function G_dn_update(model::HubbardDQMC, n::Int64, i::Int64, G::Matrix{Float64})
    a = G[:, i]
    b = ((I - G)[i, :])'
    G - Δ_dn(model, n, i) * a * b / accept_ratio_down(model, n, i, G)
end

function sweep!(model::HubbardDQMC, n_sweep::Int64)
    n_imtimes = model.n_imtimes
    s = model.s
    lattice = model.lattice
    n_sites = lattice.n_sites

    current_G_up = G_up_stable(model, 1)
    current_G_dn = G_dn_stable(model, 1)

    function sweep(τ)
        for i in 1:n_sites

            accept_rate_up = accept_ratio_up(model, τ, i, current_G_up)
            accept_rate_down = accept_ratio_down(model, τ, i, current_G_dn)
            accept_rate = accept_rate_up * accept_rate_down
            
            if rand(Float64) < accept_rate
                current_G_up = G_up_update(model, τ, i, current_G_up)
                current_G_dn = G_dn_update(model, τ, i, current_G_dn)
                s[τ, i] *= -1
            end
        end
    end

    for _ in 1:n_sweep
        for τ in 1:(n_imtimes-1)
            sweep(τ)
            
            if τ % n_wrap == 0
                current_G_up = G_dn_stable(model, τ + 1)
                current_G_dn = G_dn_stable(model, τ + 1)
            else
                current_G_up = B_up_τ(model, τ + 1) * current_G_up * inv(B_up_τ(model, τ + 1))
                current_G_dn = B_dn_τ(model, τ + 1) * current_G_dn * inv(B_dn_τ(model, τ + 1))
            end

        end

        for τ in n_imtimes:-1:2
            sweep(τ)

            if τ % n_wrap == 0
                current_G_up = G_up_stable(model, τ - 1)
                current_G_dn = G_dn_stable(model, τ - 1)
            else
                current_G_up = inv(B_up_τ(model, τ)) * current_G_up * B_up_τ(model, τ)
                current_G_dn = inv(B_dn_τ(model, τ)) * current_G_dn * B_dn_τ(model, τ)
            end

        end
    end
end