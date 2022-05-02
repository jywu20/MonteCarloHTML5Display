# Single-τ B-matrix
# Possible optimization:
# - Checkboard decomposition

function B_up(τ)
    diagm(exp.(α * s_τ[τ, :])) * exp(- Δτ * T_kin)
end

function B_up_inv(τ)
    exp(Δτ * T_kin) * diagm(- exp.(α * s_τ[τ, :]))
end

function B_dn(τ)
    diagm(exp.(- α * s_τ[τ, :])) * exp(- Δτ * T_kin)
end

function B_dn_inv(τ)
    exp(Δτ * T_kin) * diagm(exp.(α * s_τ[τ, :]))
end

function B_up_0_τ(τ)
    B = I
    for τp in 1 : τ
        B = B_up(τp) * B
    end
    B
end

function B_dn_0_τ(τ)
    B = I
    for τp in 1 : τ
        B = B_dn(τp) * B
    end
    B
end

function B_up_τ_β(τ)
    B = I
    for τp in τ+1 : n_τ
        B = B_up(τp) * B
    end
    B
end

function B_dn_τ_β(τ)
    B = I
    for τp in τ+1:n_τ
        B = B_dn(τp) * B
    end
    B
end

function udv_decompose(m::Matrix{Float64})::Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}
    U, S, Vt = svd(m)
    (U, diagm(S), Vt')
end

function B_up_0_τ_udv(τ)
    if τ == 0
        return (I, I, I)
    end
    U, D, V = udv_decompose(B_up_storage[:, :, 1])
    for τp in 2 : τ 
        Up, Dp, Vp = udv_decompose(B_up_storage[:, :, τp] * U * D)
        copy!(V, Vp * V)
        copy!(D, Dp)
        copy!(U, Up)
    end
    (U, D, V)
end

function B_dn_0_τ_udv(τ)
    if τ == 0
        return (I, I, I)
    end
    U, D, V = udv_decompose(B_down_storage[:, :, 1])
    for τp in 2 : τ
        Up, Dp, Vp = udv_decompose(B_down_storage[:, :, τp] * U * D)
        copy!(V, Vp * V)
        copy!(D, Dp)
        copy!(U, Up)
    end
    (U, D, V)
end

function B_up_τ_β_vdu(τ)
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

function B_dn_τ_β_vdu(τ)
    if τ == n_τ
        return (I, I, I)
    end
    V, D, U = udv_decompose(B_down_storage[:, :, n_τ])
    for τp in n_τ - 1 : -1 : τ + 1
        Vp, Dp, Up = udv_decompose(D * U * B_down_storage[:, :, τp])
        copy!(V, V * Vp)
        copy!(D, Dp)
        copy!(U, Up)
    end
    (V, D, U)
end