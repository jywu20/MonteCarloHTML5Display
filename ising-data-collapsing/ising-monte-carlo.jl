using Random

#region The abstract, lattice independent Ising model

"""
Generic parameters:

- `S` The type of site coordinate
- `F` The type of float number used

Fields:

- `σ` The spin configuration.
- `site_list[i]` gives the coordinates of `i`
- `inverse_list[coord]` gives the index of site with coordinates `coord` in `site_list`
- `neighbor_list[i][j]` gives the `j`th neighbor of site indexed as `i` in `site_list`
- `J` and `h` gives the following free energy
   F = - J ∑_{⟨i,j⟩} σ_i σ_j - h ∑_i σ_i
   and the partition function is 
   Z = ∑_σ exp(- β F) 
"""
struct IsingMetropolis{S, F<:AbstractFloat} 
    σ::Vector{Int}
    site_list::Vector{S}
    inverse_list::Dict{S, Int}
    neighbor_list::Vector{Vector{Int}}
    J::F
    h::F
    β::F
end

function IsingMetropolis(::Type{S}, ::Type{F}, 
    site_list::Vector{S}, inverse_list::Dict{S, Int}, neighbor_list::Vector{Vector{Int}}, 
    J::F, h::F, β::F) where {S, F}
    σ = ones(Int, length(site_list))
    if h < 0
        σ = - σ
    end
    IsingMetropolis{S, F}(σ, site_list, inverse_list, neighbor_list, J, h, β)
end

function update!(model::IsingMetropolis)
    σ = model.σ
    site_list = model.site_list
    neighbor_list = model.neighbor_list
    J = model.J
    h = model.h
    β = model.β

    for i in 1 : length(site_list)
        neighbors = neighbor_list[i]
        accept_rate_J = - 2 * β * J * sum(map(j -> σ[i] * σ[j], neighbors))
        accept_rate_h = - 2 * β * h * σ[i]
        accept_rate = exp(accept_rate_J + accept_rate_h)
        if rand() < accept_rate
            σ[i] *= -1
        end
    end
end

function magnetization(model::IsingMetropolis)
    sum(model.σ) / length(model.site_list)
end

#endregion

#region Construct the lattice

function square_lattice_2D(L::Int)::Tuple{
    Vector{Tuple{Int, Int}}, # The type of `site_list` 
    Dict{Tuple{Int, Int}, Int}} # The type of `inverse_list`
    site_list = Vector{Tuple{Int, Int}}(undef, L^2)
    inverse_list = Dict{Tuple{Int, Int}, Int}()

    site_index = 0
    for y in 1 : L
        for x in 1 : L
            site_index += 1

            site_list[site_index] = (x, y)
            push!(inverse_list, (x, y) => site_index)
        end
    end

    (site_list, inverse_list)
end

function back_into_range(period, value)
    result = (value + period) % period
    if result == 0
        return period
    end
    
    result
end

function nearest_neighbor_square_lattice_2D(L::Int)
    site_list, inverse_list = square_lattice_2D(L)
    site_num = length(site_list)
    neighbor_list = Vector{Vector{Int}}(undef, site_num)

    for i in 1 : site_num
        x_i, y_i = site_list[i]
        j1 = inverse_list[(back_into_range(L, x_i + 1), y_i)]
        j2 = inverse_list[(back_into_range(L, x_i - 1), y_i)]
        j3 = inverse_list[(x_i, back_into_range(L, y_i + 1))]
        j4 = inverse_list[(x_i, back_into_range(L, y_i - 1))]

        neighbor_list[i] = [j1, j2, j3, j4]
    end

    (site_list, inverse_list, neighbor_list)
end

#endregion