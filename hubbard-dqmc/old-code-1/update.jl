function Δ_up(τ, i)
    res = zeros(n_sites, n_sites)
    res[i, i] = exp(- 2 * α * s_τ[τ, i]) - 1
    res
end

function Δ_down(τ, i)
    res = zeros(n_sites, n_sites)
    res[i, i] = exp(2 * α * s_τ[τ, i]) - 1
    res
end

function B_up_updated(τ, i)
    (I + Δ_up(τ, i)) * B_up(τ)
end

function B_down_updated(τ, i)
    (I + Δ_down(τ, i)) * B_down(τ)
end

function accept_rate_up()
    1 + Δ_up(τ_now, i_now)[i_now, i_now] * (1 - G_up(τ_now)[i_now, i_now])
end

function accept_rate_down()
    1 + Δ_down(τ_now, i_now)[i_now, i_now] * (1 - G_down(τ_now)[i_now, i_now])
end

function accept_rate()
    accept_rate_up() * accept_rate_down()
end

function G_update()
    copy!(G_up, G_up - G_up * Δ_up(τ_now, i_now) * (I - G_up) / accept_rate_up())
    copy!(G_dn, G_dn - G_dn * Δ_down(τ_now, i_now) * (I - G_dn) / accept_rate_down())
    s_τ[τ_now, i_now] *= -1
    B_up_storage[:, :, τ_now] = B_up(τ_now)
    B_down_storage[:, :, τ_now] = B_down(τ_now)
end

function wrap()
    copy!(G_up, G_up(τ_now))
    copy!(G_dn, G_down(τ_now))
end