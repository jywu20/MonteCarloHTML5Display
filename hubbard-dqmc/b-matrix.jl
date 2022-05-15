# Single-τ B-matrix
# Possible optimization:
# - Checkboard decomposition

function B_up(model::HubbardDQMC, τ)
    α = model.α
    s_τ = model.s
    Δτ = model.Δτ
    T_kin = model.T
    diagm(exp.(α * s_τ[τ, :])) * exp(- Δτ * T_kin)
end

function B_up_inv(model::HubbardDQMC, τ)
    α = model.α
    s_τ = model.s
    Δτ = model.Δτ
    T_kin = model.T
    exp(Δτ * T_kin) * diagm(- exp.(α * s_τ[τ, :]))
end

function B_dn(model::HubbardDQMC, τ)
    α = model.α
    s_τ = model.s
    Δτ = model.Δτ
    T_kin = model.T
    diagm(exp.(- α * s_τ[τ, :])) * exp(- Δτ * T_kin)
end

function B_dn_inv(model::HubbardDQMC, τ)
    α = model.α
    s_τ = model.s
    Δτ = model.Δτ
    T_kin = model.T
    exp(Δτ * T_kin) * diagm(exp.(α * s_τ[τ, :]))
end

function B_up_0_τ(model::HubbardDQMC, τ)
    B = I 
    for τ′ in 1 : τ
        B = B_up(model, τ′) * B
    end
    B
end

function B_dn_0_τ(model::HubbardDQMC, τ)
    B = I
    for τ′ in 1 : τ
        B = B_dn(model, τ′) * B
    end
    B
end

function B_up_τ_β(model::HubbardDQMC, τ)
    B = I
    for τ′ in τ+1 : n_τ
        B = B_up(model, τ′) * B
    end
    B
end

function B_dn_τ_β(model::HubbardDQMC, τ)
    B = I
    for τ′ in τ+1:n_τ
        B = B_dn(model, τ′) * B
    end
    B
end

function udv_decompose(m::Matrix{Float64})::Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}
    U, S, Vt = svd(m)
    (U, diagm(S), Vt')
end

function B_up_0_τ_udv(model::HubbardDQMC, τ)
    B_up_storage = model.B_up_storage
    if τ == 0
        return (I, I, I)
    end
    U, D, V = udv_decompose(B_up_storage[:, :, 1])
    for τ′ in 2 : τ 
        Up, Dp, Vp = udv_decompose(B_up_storage[:, :, τ′] * U * D)
        copy!(V, Vp * V)
        copy!(D, Dp)
        copy!(U, Up)
    end
    (U, D, V)
end

function B_dn_0_τ_udv(model::HubbardDQMC, τ)
    B_dn_storage = model.B_dn_storage
    if τ == 0
        return (I, I, I)
    end
    U, D, V = udv_decompose(B_dn_storage[:, :, 1])
    for τp in 2 : τ
        Up, Dp, Vp = udv_decompose(B_dn_storage[:, :, τp] * U * D)
        copy!(V, Vp * V)
        copy!(D, Dp)
        copy!(U, Up)
    end
    (U, D, V)
end

function B_up_τ_β_vdu(model::HubbardDQMC, τ)
    B_up_storage = model.B_up_storage
    if τ == n_τ
        return (I, I, I)
    end
    V, D, U = udv_decompose(B_up_storage[:, :, n_τ])
    for τp in n_τ - 1 : -1 : τ + 1
        Vp, Dp, Up = udv_decompose(D * U * B_up_storage[:, :, τp])
        copy!(V, V * Vp)
        copy!(D, Dp)
        copy!(U, Up)
    end
    (V, D, U)
end

function B_dn_τ_β_vdu(model::HubbardDQMC, τ)
    B_dn_storage = model.B_dn_storage
    if τ == n_τ
        return (I, I, I)
    end
    V, D, U = udv_decompose(B_dn_storage[:, :, n_τ])
    for τp in n_τ - 1 : -1 : τ + 1
        Vp, Dp, Up = udv_decompose(D * U * B_dn_storage[:, :, τp])
        copy!(V, V * Vp)
        copy!(D, Dp)
        copy!(U, Up)
    end
    (V, D, U)
end