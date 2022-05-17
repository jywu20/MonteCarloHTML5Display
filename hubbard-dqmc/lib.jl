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

    lattice::L

    # s[τ, i]
    s::Matrix{Int64}

    # T[i, j]
    T::Matrix{Float64}
    # exp(Δτ * T)
    expT::Matrix{Float64}
    # exp(- Δτ * T)
    expmT::Matrix{Float64}

    id::Matrix{Float64}
    B_up_storage::Array{Float64, 3}
    B_dn_storage::Array{Float64, 3}
end

include("b-matrix.jl")
include("green.jl")
include("update.jl")
include("init.jl")
include("sweep.jl")