#region Definitions

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
    exp(Δτ * T_kin) * diagm(exp.(- α * s_τ[τ, :]))
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

#endregion

#region Multiplication of exp(±V)

"""
In-situ diagm(exp.(α * s_τ[τ, :])) * m
"""
function lmul_exp_V!(model::HubbardDQMC, τ, m)
    α = model.α
    s_τ = model.s
    lmul!(Diagonal(exp.(α * s_τ[τ, :])), m)
end

"""
In-situ diagm(exp.(- α * s_τ[τ, :])) * m
"""
function lmul_exp_mV!(model::HubbardDQMC, τ, m)
    α = model.α
    s_τ = model.s
    lmul!(Diagonal(exp.(- α * s_τ[τ, :])), m)
end

"""
In-situ m * diagm(exp.(α * s_τ[τ, :]))
"""
function rmul_exp_V!(model::HubbardDQMC, τ, m)
    α = model.α
    s_τ = model.s
    rmul!(m, Diagonal(exp.(α * s_τ[τ, :])))
end

"""
In-situ m * diagm(exp.(- α * s_τ[τ, :]))
"""
function rmul_exp_mV!(model::HubbardDQMC, τ, m)
    α = model.α
    s_τ = model.s
    rmul!(m, Diagonal(exp.(- α * s_τ[τ, :])))
end

#endregion

#region Multiplication of exp(±T)

"""
In-situ exp(Δτ * T) * m
"""
function lmul_exp_T!(model::HubbardDQMC, m)
    copy!(m, model.expT * m)
end

"""
In-situ exp(- Δτ * T) * m
"""
function lmul_exp_mT!(model::HubbardDQMC, m)
    copy!(m, model.expmT * m)
end

"""
In-situ m * exp(Δτ * T) 
"""
function rmul_exp_T!(model::HubbardDQMC, m)
    copy!(m, m * model.expT)
end

"""
In-situ m * exp(- Δτ * T) 
"""
function rmul_exp_mT!(model::HubbardDQMC, m)
    copy!(m, m * model.expmT)
end

#endregion

#region 

function lmul_B_up!(model::HubbardDQMC, τ, m)
    lmul_exp_mT!(model, m)
    lmul_exp_V!(model, τ, m)
end

function lmul_B_dn!(model::HubbardDQMC, τ, m)
    lmul_exp_mT!(model, m)
    lmul_exp_mV!(model, τ, m)
end

function lmul_B_up_inv!(model::HubbardDQMC, τ, m)
    lmul_exp_mV!(model, τ, m)
    lmul_exp_T!(model, m)
end

function lmul_B_dn_inv!(model::HubbardDQMC, τ, m)
    lmul_exp_V!(model, τ, m)
    lmul_exp_T!(model, m)
end

function rmul_B_up!(model::HubbardDQMC, τ, m)
    rmul_exp_V!(model, τ, m)
    rmul_exp_mT!(model, m)
end

function rmul_B_dn!(model::HubbardDQMC, τ, m)
    rmul_exp_mV!(model, τ, m)
    rmul_exp_mT!(model, m)
end

function rmul_B_up_inv!(model::HubbardDQMC, τ, m)
    rmul_exp_T!(model, m)
    rmul_exp_mV!(model, τ, m)
end

function rmul_B_dn_inv!(model::HubbardDQMC, τ, m)
    rmul_exp_T!(model, m)
    rmul_exp_V!(model, τ, m)
end

#endreion