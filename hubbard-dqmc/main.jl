# The main Monte Carlo loop

include("lib.jl")

L = 4 
U = 4.0
β = 4.0
n_imtimes = 100
n_wrap = 10

lattice = SquareLattice2D(L)
model = HubbardDQMC(SquareLattice2D, lattice, U, β, n_imtimes, n_wrap)
sweep!(model, 10)