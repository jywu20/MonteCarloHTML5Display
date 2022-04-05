using Statistics
include("ising-monte-carlo.jl")

Ls = [10, 20, 30, 40, 50]
ΔT = 0.05
T_min = 2.0
T_max = 2.4
Ts = T_min : ΔT : T_max

heating_up_time = 2000
recording_time = 5000
magnetization_record = zeros(recording_time)

working_path = "D:\\Projects\\Modern Physics Experiments\\HTML5\\ising-data-collapsing\\"
name = "finite-size-scaling-data-1.csv"

open(working_path * name, "w") do file
    println(file, "L, T, m, χ")
end

J = 1.0
h = 0.0

for L in Ls
    site_list, inverse_list, neighbor_list = nearest_neighbor_square_lattice_2D(L)
    for T in Ts
        β = 1 / T
        model = IsingMetropolis(Tuple{Int, Int}, Float64, site_list, inverse_list, neighbor_list, J, h, β)

        for _ in 1 : heating_up_time
            update!(model)
        end

        for update_count in 1 : recording_time
            update!(model)
            magnetization_record[update_count] = abs(magnetization(model))
        end

        m = mean(magnetization_record)
        m² = mean(magnetization_record.^2)
        χ = β * L^2 * (m² - m^2)

        open(working_path * name, "a") do file
            println(file, "$L, $T, $m, $χ")
        end
    end
end