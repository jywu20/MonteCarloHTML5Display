# Calculate some intermediate contstants.

# The time step in path integral
Δτ = β / n_τ

α = acosh(exp(Δτ * U / 2))

# We are working on a two dimensional lattice.
n_sites = n_side^2

log_path = path * "log"
observables_path = path * "observables"
error_path = path * "error"

#####################################################################################################
# The kinetic Hamiltonian
T_kin = zeros((n_sites, n_sites))
for i in 1:n_sites
    for j in nearest_neighbors(i)
        T_kin[i, j] = 1.0
    end
end
T_kin = - t * T_kin

#####################################################################################################
# Initialize the system configuration.

# What imaginary time step are we at now. Ranging from 1 to n_τ
τ_now = 1
# Which site are we on now. Ranging from 1 to n_sites
i_now = 1

# The auxiliary field. The first index is the imaginary time, the second index labels sites.
s_τ = rand((-1, 1), n_τ, n_sites)

# Stroed B-matrices, to avoid unncessisary calculations.
B_up_storage = zeros(n_sites, n_sites, n_τ)
B_down_storage = zeros(n_sites, n_sites, n_τ)

for τ in 1 : n_τ
    B_up_storage[:, :, τ] = B_up(τ)
    B_down_storage[:, :, τ] = B_down(τ)
end

# The Green functions of the current time τ.
G_up = G_up(τ_now)
G_dn = G_down(τ_now)

# Counts how many times have `propag_forward` and `propag_backward` been invoked.
# If it reaches n_wrap, then Green functions will be calculated from B-matrices, and the counter is set back to 0.
wrap_count = 0

#####################################################################################################
# Logging.

function init_logging()
    open(log_path, "w") do log_file
        println(log_file, "Welcome to a simple DQMC simulation of Hubbard model.")
        println(log_file)

        # Print parameters used in the simulation.
        println(log_file, "We are simulating with")
        println(log_file, "  t      =   $(string(t))")
        println(log_file, "  U      =   $(string(U))")
        println(log_file, "  beta   =   $(string(β))")
        println(log_file, "  dtau   =   $(string(Δτ))")
        println(log_file, "  alpha  =   $(string(α))")
        println(log_file, "on a $n_side * $n_side 2d lattice and with $n_τ time steps.")

        if double_check
            println(log_file)
            println(log_file, "Double checking is on. During the calculation, the program double-checks the error between the optimized version of B-matrices, Green functions, etc. and their by-definition values.")
            println(log_file, "This can be extremly slow. Turn off the marker when scaling up.")
            println(log_file)
        end

        println(log_file, "Initialization completed.")
    end
    
    if double_check
        open(error_path, "w") do error_file
            println(error_file, "numerical error found in the double-check process")
        end
    end

    open(observables_path, "w") do observable_file
        println(observable_file, "Observables:")
    end
end

function dqmc_log(content::String)
    open(log_path, "a") do log_file
        println(log_file, content)
    end
end