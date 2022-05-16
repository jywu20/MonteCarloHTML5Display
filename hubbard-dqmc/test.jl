using Plots

include("lib.jl")
include("config.jl")

##

lattice = SquareLattice2D(n_side)

n_sites = lattice.n_sites
s = rand((-1, 1), (n_τ, n_sites))

model = HubbardDQMC(SquareLattice2D, lattice, t, U, β, n_τ, n_wrap, s = s)

##

α = model.α
s_τ = model.s
τ = 4
id = model.id
A = rand(size(id)...)

@time diagm(exp.(- α * s_τ[τ, :])) * A 
##

@time lmul_exp_mV!(model, τ, A)

## 
