function G_up(τ)
    Ur, Dr, Vr = B_up_0_τ_udv(τ)
    Vl, Dl, Ul = B_up_τ_β_vdu(τ)
    U, D, V = udv_decompose(inv(Ul * Ur) + Dr * (Vr * Vl) * Dl)
    inv(V * Ul) * inv(D) * inv(Ur * U)
end

function G_down(τ)
    Ur, Dr, Vr = B_down_0_τ_udv(τ)
    Vl, Dl, Ul = B_down_τ_β_vdu(τ)
    U, D, V = udv_decompose(inv(Ul * Ur) + Dr * (Vr * Vl) * Dl)
    inv(V * Ul) * inv(D) * inv(Ur * U)
end

function propag_forward()
    global wrap_count += 1
    global τ_now += 1
    copy!(G_up, B_up(τ_now) * G_up * inv(B_up(τ_now)))
    copy!(G_dn, B_down(τ_now) * G_dn * inv(B_down(τ_now)))

    if wrap_count == n_wrap
        global wrap_count = 0
        copy!(G_up, G_up(τ_now))
        copy!(G_dn, G_down(τ_now))
    end
end

function propag_backward()
    global wrap_count += 1
    copy!(G_up, B_up_inv(τ_now) * G_up * B_up(τ_now))
    copy!(G_dn, B_down_inv(τ_now) * G_dn * B_down(τ_now))
    global τ_now -= 1

    if wrap_count == n_wrap
        global wrap_count = 0
        copy!(G_up, G_up(τ_now))
        copy!(G_dn, G_down(τ_now))
        return
    end
end