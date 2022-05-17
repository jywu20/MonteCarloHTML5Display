# This file contains constants that determines the model Hamiltonian and controls the simulation process.

# Constants in model Hamiltonian
# The model Hamiltonian is H = - t ∑_{i, σ} c†_{iσ} c_{iσ} + U \sum_i (n_{i↑} - 1/2) (n_{i↓} - 1/2)
t = 1.0
U = 4.0

#####################################################################################################

# Constants controlling the simulation

# The inverse temperature
β = 4.0
# The number of time steps in path integral
n_τ = 40

# The side length of the lattice
n_side = 4

# Controlling heating-up times, sweeping times and bin number
heating_up_steps = 0
n_bin = 1
n_sweep = 100
n_wrap = 20

# Decide whether error tracking is performed
double_check = false 

#####################################################################################################

# Whether to show progress bar
show_progress = true

# Whether to only show binned data.
show_binned_only = false

# Whether to show standard error in each bin.
show_bin_std = true 

#####################################################################################################
