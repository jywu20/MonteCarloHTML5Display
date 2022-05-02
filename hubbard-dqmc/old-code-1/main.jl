using LinearAlgebra
using ProgressMeter
using Statistics

include("config.jl")
include("lattice.jl")
include("b-matrix.jl")
include("green.jl")
include("update.jl")

include("observe.jl")

include("init.jl")

if init_log
    init_logging()
end

include("core.jl")

dqmc_log("Heating up for $heating_up_steps steps.")
sweep(heating_up_steps, true)
dqmc_log("Heating up completed.")
dqmc_log("")

for bin_count in 1 : n_bin
    sweep(n_sweep, false)
    binning()
    dqmc_log("Bin $bin_count completed.")
end