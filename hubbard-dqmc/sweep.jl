function sweep!(model::HubbardDQMC, n_sweep::Int64; observe = nothing)
    # Counts how many times have `propag_forward` and `propag_backward` been invoked.
    # If it reaches n_wrap, then Green functions will be calculated from B-matrices, and the counter is set back to 0.
    wrap_count = 0

    # What imaginary time step are we at now. Ranging from 1 to n_τ
    τ_now = 1

    # The Green functions of the current time τ.
    G_up = G_up_τ(model, τ_now)
    G_dn = G_dn_τ(model, τ_now)

    n_sites = model.lattice.n_sites

    for _ in 1 : n_sweep
        for _ in 1 : n_τ - 1
            for i in 1 : n_sites
                if rand() < accept_rate(model, G_up, G_dn, τ_now, i)
                    G_update(model, G_up, G_dn, τ_now, i)
                end
            end
            
            # Propagate forward.
            wrap_count += 1
            τ_now += 1
            #copy!(G_up, B_up(model, τ_now) * G_up * B_up_inv(model, τ_now))
            lmul_B_up!(model, τ_now, G_up)
            rmul_B_up_inv!(model, τ_now, G_up)
            #copy!(G_dn, B_dn(model, τ_now) * G_dn * B_dn_inv(model, τ_now))
            lmul_B_dn!(model, τ_now, G_dn)
            rmul_B_dn_inv!(model, τ_now, G_dn)

            if wrap_count == n_wrap
                wrap_count = 0
                copy!(G_up, G_up_τ(model, τ_now))
                copy!(G_dn, G_dn_τ(model, τ_now))
            end
        end

        if observe !== nothing
            observe(model, G_up, G_dn)
        end
        
        for _ in 1 : n_τ - 1
            for i in 1 : n_sites
                if rand() < accept_rate(model, G_up, G_dn, τ_now, i)
                    G_update(model, G_up, G_dn, τ_now, i)
                end
            end
            
            # Propagate backward
            wrap_count += 1
            # copy!(G_up, B_up_inv(model, τ_now) * G_up * B_up(model, τ_now))
            lmul_B_up_inv!(model, τ_now, G_up)
            rmul_B_up!(model, τ_now, G_up)
            # copy!(G_dn, B_dn_inv(model, τ_now) * G_dn * B_dn(model, τ_now))
            lmul_B_dn_inv!(model, τ_now, G_dn)
            rmul_B_dn!(model, τ_now, G_dn)
            τ_now -= 1
        
            if wrap_count == n_wrap
                wrap_count = 0
                copy!(G_up, G_up_τ(model, τ_now))
                copy!(G_dn, G_dn_τ(model, τ_now))
            end
        end
    end
end
