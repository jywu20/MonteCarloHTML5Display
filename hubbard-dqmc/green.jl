function G_up_τ(model::HubbardDQMC, τ)
    Ur, Dr, Vr = B_up_0_τ_udv(model, τ)
    Vl, Dl, Ul = B_up_τ_β_vdu(model, τ)
    U, D, V = udv_decompose(inv(Ul * Ur) + Dr * (Vr * Vl) * Dl)
    inv(V * Ul) * inv(D) * inv(Ur * U)
end

function G_dn_τ(model::HubbardDQMC, τ)
    Ur, Dr, Vr = B_dn_0_τ_udv(model, τ)
    Vl, Dl, Ul = B_dn_τ_β_vdu(model, τ)
    U, D, V = udv_decompose(inv(Ul * Ur) + Dr * (Vr * Vl) * Dl)
    inv(V * Ul) * inv(D) * inv(Ur * U)
end
