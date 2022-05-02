function Δ_up(τ, i)
    res = zeros(n_sites, n_sites)
    res[i, i] = exp(- 2 * α * s_τ[τ, i]) - 1
    res
end

function Δ_dn(τ, i)
    res = zeros(n_sites, n_sites)
    res[i, i] = exp(2 * α * s_τ[τ, i]) - 1
    res
end

function B_up_updated(τ, i)
    (I + Δ_up(τ, i)) * B_up(τ)
end

function B_dn_updated(τ, i)
    (I + Δ_dn(τ, i)) * B_dn(τ)
end

function accept_rate_up(τ, i)
    1 + Δ_up(τ, i)[i, i] * (1 - G_up_τ(τ)[i, i])
end

function accept_rate_down(τ, i)
    1 + Δ_dn(τ, i)[i, i] * (1 - G_dn_τ(τ)[i, i])
end

function accept_rate(τ, i)
    accept_rate_up(τ, i) * accept_rate_down(τ, i)
end

function G_update(τ, i)
    copy!(G_up, G_up - G_up * Δ_up(τ, i) * (I - G_up) / accept_rate_up(τ, i))
    copy!(G_dn, G_dn - G_dn * Δ_dn(τ, i) * (I - G_dn) / accept_rate_down(τ, i))
    s_τ[τ, i] *= -1
    B_up_storage[:, :, τ] = B_up(τ)
    B_down_storage[:, :, τ] = B_dn(τ)
end

function wrap()
    copy!(G_up, G_up_τ(τ_now))
    copy!(G_dn, G_dn_τ(τ_now))
end