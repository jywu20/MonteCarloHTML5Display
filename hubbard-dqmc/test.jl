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

@time sweep_old!(model, 4)
@time sweep!(model, 4)
nothing

##

let 
    id = copy(model.id)
    τ = 7
    A1 = rand(size(id)...)
    A2 = copy(A1)
   
    @time begin
        lmul_B_up!(model, τ, A1)
        rmul_B_up_inv!(model, τ, A1)
    end 
    @time B_up(model, τ) * A2 * B_up_inv(model, τ)
    nothing
end

##

@profview sweep!(model, 4)

##
@time sweep!(model, 4)

