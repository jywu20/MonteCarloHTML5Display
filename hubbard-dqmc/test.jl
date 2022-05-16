using Plots

include("lib.jl")
include("config.jl")

##

lattice = SquareLattice2D(n_side)

n_sites = lattice.n_sites
s = rand((-1, 1), (n_τ, n_sites))

model = HubbardDQMC(SquareLattice2D, lattice, t, U, β, n_τ, n_wrap, s = s)

##

"""
In-situ diagm(exp.(α * s_τ[τ, :])) * m
"""
function lmul_exp_V!(model::HubbardDQMC, τ, m)
    α = model.α
    s_τ = model.s
    lmul!(Diagonal(exp.(α * s_τ[τ, :])), m)
end

"""
In-situ m * diagm(exp.(- α * s_τ[τ, :]))
"""
function rmul_exp_mV!(model::HubbardDQMC, τ, m)
    α = model.α
    s_τ = model.s
    rmul!(m, Diagonal(exp.(- α * s_τ[τ, :])))
end



## 
# The speed of diagonal product 
