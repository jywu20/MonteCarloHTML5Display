# Calculating B-matrices by definition.

#####################################################################################################

# Single-τ B-matrix
# Possible optimization:
# - Checkboard decomposition

function B_up_def(τ)
    diagm(exp.(α * s_τ[τ, :])) * exp(- Δτ * T_kin)
end

function B_down_def(τ)
    diagm(exp.(- α * s_τ[τ, :])) * exp(- Δτ * T_kin)
end

#####################################################################################################

# Double time B-matrix
# Possible optimization:
# - Wrapping, or numerical stablization

function B_up_0_τ_def(τ)
    B = I
    for τp in 1:τ
        B = B_up_def(τp) * B
    end
    B
end

function B_down_0_τ_def(τ)
    B = I
    for τp in 1:τ
        B = B_down_def(τp) * B
    end
    B
end

function B_up_τ_β_def(τ)
    B = I
    for τp in τ+1:n_τ
        B = B_up_def(τp) * B
    end
    B
end

function B_down_τ_β_def(τ)
    B = I
    for τp in τ+1:n_τ
        B = B_down_def(τp) * B
    end
    B
end

function G_up_def(τ)
    inv(I + B_up_0_τ_def(τ) * B_up_τ_β_def(τ))
end

function G_down_def(τ)
    inv(I + B_down_0_τ_def(τ) * B_down_τ_β_def(τ))
end