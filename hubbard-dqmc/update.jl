function Δ_up(model::HubbardDQMC, τ, i)
    α = model.α
    s_τ = model.s
    n_sites = model.lattice.n_sites
    res = zeros(n_sites, n_sites)
    res[i, i] = exp(- 2 * α * s_τ[τ, i]) - 1
    res
end

function Δ_dn(model::HubbardDQMC, τ, i)
    α = model.α
    s_τ = model.s
    n_sites = model.lattice.n_sites
    res = zeros(n_sites, n_sites)
    res[i, i] = exp(2 * α * s_τ[τ, i]) - 1
    res
end

function B_up_updated(model::HubbardDQMC, τ, i)
    (I + Δ_up(model, τ, i)) * B_up(model, τ)
end

function B_dn_updated(model::HubbardDQMC, τ, i)
    (I + Δ_dn(model, τ, i)) * B_dn(model, τ)
end

function accept_rate_up(model::HubbardDQMC, τ, i)
    1 + Δ_up(model, τ, i)[i, i] * (1 - G_up_τ(model, τ)[i, i])
end

function accept_rate_dn(model::HubbardDQMC, τ, i)
    1 + Δ_dn(model, τ, i)[i, i] * (1 - G_dn_τ(model, τ)[i, i])
end

function accept_rate(model::HubbardDQMC, τ, i)
    accept_rate_up(model, τ, i) * accept_rate_dn(model, τ, i)
end

function G_update(model::HubbardDQMC, G_up, G_dn, τ, i)
    s_τ = model.s
    B_up_storage = model.B_up_storage
    B_dn_storage = model.B_dn_storage
    copy!(G_up, G_up - G_up * Δ_up(model, τ, i) * (I - G_up) / accept_rate_up(model, τ, i))
    copy!(G_dn, G_dn - G_dn * Δ_dn(model, τ, i) * (I - G_dn) / accept_rate_dn(model, τ, i))
    s_τ[τ, i] *= -1
    B_up_storage[:, :, τ] = B_up(model, τ)
    B_dn_storage[:, :, τ] = B_dn(model, τ)
end

function wrap()
    copy!(G_up, G_up_τ(model, τ_now))
    copy!(G_dn, G_dn_τ(model, τ_now))
end