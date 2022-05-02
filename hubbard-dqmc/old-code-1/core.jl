function sweep(times::Int64, heating::Bool)
    progress = Progress(times)

    for _ in 1 : times
        for _ in 1 : n_τ - 1
            for i in 1 : n_sites
                global i_now = i
                if rand() < accept_rate()
                    G_update()
                end
            end
            propag_forward()    
        end

        if ! heating
            observe()
        end
        
        for _ in 1 : n_τ - 1
            for i in 1 : n_sites
                global i_now = i
                if rand() < accept_rate()
                    G_update()
                end
            end
            propag_backward()    
        end
    
        if show_progress
            next!(progress)
        end
    end
end
