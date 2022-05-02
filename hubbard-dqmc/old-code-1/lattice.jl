function site_to_coord(i::Int64)::Tuple{Int64, Int64}
    y = i % n_side
    if y == 0
        y = n_side
    end
    x = 1 + Int((i - y) / n_side)
    
    (x, y)
end

function coord_to_site(x::Int64, y::Int64)::Int64
    (x - 1) * n_side + y
end

function nearest_neighbors(i::Int64)::Vector{Int64}
    x, y = site_to_coord(i)
    [
        coord_to_site((i -> i > 0 ? i : i + n_side)(mod(x + 1, n_side)), y),
        coord_to_site((i -> i > 0 ? i : i + n_side)(mod(x - 1, n_side)), y),
        coord_to_site(x, (i -> i > 0 ? i : i + n_side)(mod(y + 1, n_side))),
        coord_to_site(x, (i -> i > 0 ? i : i + n_side)(mod(y - 1, n_side)))
    ]
end