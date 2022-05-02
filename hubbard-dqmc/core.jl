function sweep(n_sweep::Int64, heating::Bool)
    progress = Progress(n_sweep)

    for _ in 1 : n_sweep
        for _ in 1 : n_τ - 1
            for i in 1 : n_sites
                if rand() < accept_rate(τ_now, i)
                    G_update(τ_now, i)
                end
            end
            propag_forward()    
        end

        if ! heating
            observe()
        end
        
        for _ in 1 : n_τ - 1
            for i in 1 : n_sites
                if rand() < accept_rate(τ_now, i)
                    G_update(τ_now, i)
                end
            end
            propag_backward()    
        end
    
        if show_progress
            next!(progress)
        end
    end
end
