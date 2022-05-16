using LinearAlgebra

include("utils.jl")
include("lattice.jl")

struct HubbardDQMC{L <: AbstractLattice} 
    t::Float64
    U::Float64
    β::Float64
    Δτ::Float64
    n_imtimes::Int64
    n_wrap::Int64
    α::Float64
    # Currently we don't put the Green functions in the model struct, because in the sweeping process,
    # they are dependent on τ, which shouldn't be placed into the model struct, or otherwise the code 
    # will be rather tedious
    # G[i, j, τ]
    #G_up::Array{Float64, 3}
    #G_dn::Array{Float64, 3}
    lattice::L
    # s[τ, x]
    s::Matrix{Int64}
    # T[i, j]
    T::Matrix{Float64}
    id::Matrix{Float64}
    B_up_storage::Array{Float64, 3}
    B_dn_storage::Array{Float64, 3}
end

include("b-matrix.jl")
include("green.jl")
include("update.jl")
include("init.jl")
include("sweep.jl")