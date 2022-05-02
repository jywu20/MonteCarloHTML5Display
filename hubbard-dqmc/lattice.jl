# Constructing lattices

"""
A lattice is defined by a mapping from site index to site coordinates (`site_list`) and 
its inverse (`inverse_list`), as well as its nearest neighbor list, next nearest neighbor list, etc.
We assume all abstract lattice has following fields:
- `site_list`, `site_list[i, :]` gives the coordinates of the site with index `i` (henceforth "site `i`")
- `inverse_list`, `inverse_list[site_list[i, :]...] == i`. Another choice is that we use `Tuple`s to record 
   the coordinates, because it's possible that we work on triangular lattices, so the possible values of 
   `site_list[i]` can't form a square, so if we declare `inverse_list` as `inverse_list[site_list[i]...] == i`, this may cause waste of space. We use the current convention just to keep consistent with existing codes.
- `neighbor_list`, `neighbor_list[i, n]` gives the n-th neighbor site of the site with index `i`. It's possible 
   that both nearest neighbors and next nearest neighbors are included. In this case, an optional field 
   `neighbor_list_indices` can be used to store the distance between the n-th neighbor site of the site with 
   index `i` and site `i` itself.
- `neighbor_list_indices`, if `n` is in `neighbor_list_indices[1]`, then `neighbor_list[i, n]` gives one of the 
    nearest neighbors of site `i`. If `n` is in `neighbor_list_indices[2]`, then `neighbor_list[i, n]` gives one 
    of the next nearest neighbors of site `i`, etc.
"""
abstract type AbstractLattice end

const square_lattice_2D_neighbor_list_indices = [[1, 2, 3, 4], [5, 6, 7, 8]]

struct SquareLattice2D <: AbstractLattice
    n_sites::Int

    # Two dimensional, site_list[i, 1] = x, site_list[i, 2] = y
    site_list::Matrix{Int}

    # inverse_list[x, y] = i
    inverse_list::Matrix{Int}

    # neighbor_list[i, 1:4] gives the nearest neighbors, while 
    neighbor_list::Matrix{Int}

    neighbor_list_indices
end

function SquareLattice2D(L::Integer)
    n_sites = L^2
    inverse_list = zeros(Int64, (L, L))
    site_list = zeros(Int64, (n_sites, 2))

    for i in 1:L
        for j in 1:L
            site_list[(i-1)*L+j, :] = [i, j]
            inverse_list[i, j] = (i - 1) * L + j
        end
    end

    neighbor_list = zeros(Int64, (n_sites, 8))
    for i in 1:L
        for j in 1:L
            neighbor_list[inverse_list[i, j], 1] = inverse_list[i, back_into_range(j+1, L)]
            neighbor_list[inverse_list[i, j], 2] = inverse_list[back_into_range(i+1, L), j]
            neighbor_list[inverse_list[i, j], 3] = inverse_list[i, back_into_range(j-1, L)]
            neighbor_list[inverse_list[i, j], 4] = inverse_list[back_into_range(i-1, L), j]
            neighbor_list[inverse_list[i, j], 5] = inverse_list[back_into_range(i+1, L), back_into_range(j+1, L)]
            neighbor_list[inverse_list[i, j], 6] = inverse_list[back_into_range(i+1, L), back_into_range(j-1, L)]
            neighbor_list[inverse_list[i, j], 7] = inverse_list[back_into_range(i-1, L), back_into_range(j-1, L)]
            neighbor_list[inverse_list[i, j], 8] = inverse_list[back_into_range(i-1, L), back_into_range(j+1, L)]
        end
    end

    SquareLattice2D(n_sites, site_list, inverse_list, neighbor_list, square_lattice_2D_neighbor_list_indices)
end
