function back_into_range(idx, upper)
    if idx > upper
        return idx % upper
    end
    (idx - upper) % upper + upper
end

function exp_two_hot_sym_mat(val, pos, id)
    result = copy(id)
    i, j = pos
    result[i, i] = cosh(val)
    result[j, j] = cosh(val)
    result[i, j] = sinh(val)
    result[j, i] = sinh(val)
    result
end