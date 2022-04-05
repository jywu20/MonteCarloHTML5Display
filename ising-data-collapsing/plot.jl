using PGFPlots, LaTeXStrings

##
# Data from finite-size-scaling-data-1.csv
Ts = [2, 2.01, 2.02, 2.03, 2.04, 2.05, 2.06, 2.07, 2.08, 2.09, 2.1, 2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.2, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.3, 2.31, 2.32, 2.33, 2.34, 2.35, 2.36, 2.37, 2.38, 2.39, 2.4]

χ_L_50 = [0.37 ,0.4 ,0.41 ,0.54 ,0.51 ,0.55 ,0.6 ,0.67 ,0.71 ,0.93 ,0.82 ,1.01 ,1 ,1.49 ,1.89 ,1.92 ,1.94 ,2.37 ,2.24 ,2.97 ,8.2 ,4.92 ,7.65 ,26.44 ,11.4 ,10.41 ,48.15 ,22.08 ,31.7 ,39.02 ,34.5 ,49.7 ,37.72 ,42.18 ,34.98 ,35.74 ,36.15 ,23.82 ,23.81 ,28.58 ,19.95]

working_path = "D:\\Projects\\Modern Physics Experiments\\HTML5\\ising-data-collapsing\\"
name = "finite-size-scaling-L-50-dT-0.01-run-1.pdf"

p = Axis([
        Plots.Linear(Ts, χ_L_50, onlyMarks=true),
    ],
    xlabel = L"T",
    ylabel = L"\chi",
    xmin = min(Ts...), xmax = max(Ts...), ymin = 0,
    style = "axis background/.style={fill=white}"
)

save(working_path * name, p)

##

Ts = [
    2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4
]

χ_L_50 = [
    0.39, 0.62, 1.05, 1.42, 4.44, 12.84, 51.96, 38.48, 21.97
]

working_path = "D:\\Projects\\Modern Physics Experiments\\HTML5\\ising-data-collapsing\\"
name = "finite-size-scaling-L-50-dT-0.05-run-1.pdf"

p = Axis([
        Plots.Linear(Ts, χ_L_50, onlyMarks=true),
    ],
    xlabel = L"T",
    ylabel = L"\chi",
    xmin = min(Ts...), xmax = max(Ts...), ymin = 0,
    style = "axis background/.style={fill=white}"
)

save(working_path * name, p)

##

T_c = 2.2691853
γ = 7 / 4
ν = 1

Ts = [
    2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4
]

χ_L_50 = [
    0.39, 0.62, 1.05, 1.42, 4.44, 12.84, 51.96, 38.48, 21.97
]

χ_L_40 = [
    0.41 , 0.48 , 0.85 , 1.25 , 3.43 , 6.64 , 23.12 , 24.46 , 18.08
]

χ_L_30 = [
    0.38, 0.52, 0.76, 2.22, 2.74, 10.58, 19.32, 19.04, 16.05
]

χ_L_20 = [
    0.41, 0.54, 1.09, 1.98, 2.21, 4.84, 7.87, 8.16, 8.46
]

ts = (Ts .- T_c) / T_c

p = Axis([
        Plots.Linear(ts * 20^(1 / ν), χ_L_20 * 20^(- γ / ν), onlyMarks = true, mark = "o", legendentry = L"L=20"),
        Plots.Linear(ts * 30^(1 / ν), χ_L_30 * 30^(- γ / ν), onlyMarks = true, mark = "o", legendentry = L"L=30"),
        Plots.Linear(ts * 40^(1 / ν), χ_L_40 * 40^(- γ / ν), onlyMarks = true, mark = "o", legendentry = L"L=40"),
        Plots.Linear(ts * 50^(1 / ν), χ_L_50 * 50^(- γ / ν), onlyMarks = true, mark = "o", legendentry = L"L=50"),
    ],
    xlabel = L"t L^{1 / \nu}",
    ylabel = L"\chi L^{- \gamma / \nu}",
    style = "axis background/.style={fill=white}",
    legendStyle = "draw=none",
    legendPos="north west"
)

working_path = "D:\\Projects\\Modern Physics Experiments\\HTML5\\ising-data-collapsing\\"
name = "ising-data-collapsing-run-1.pdf"
save(working_path * name, p)

##

T_c = 2.2691853
γ = 7 / 4
ν = 1

Ts = [
    2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4
]

χ_L_20 = [
    0.47, 1, 1.21, 1.08, 3.2, 4.48, 9.23, 10.25, 8.5
]

χ_L_30 = [
    0.42, 0.51, 0.95, 1.81, 2.83, 8.41, 16.4, 16.16, 14.46
]

χ_L_40 = [
    0.35, 0.52, 1.09, 2.16, 3.46, 12.93, 25.67, 30.62, 19.36
]

χ_L_50 = [
    0.37, 0.52, 0.83, 1.69, 3.01, 41.69, 40.51, 36.11, 17.07
]

ts = (Ts .- T_c) / T_c

p = Axis([
        Plots.Linear(ts * 20^(1 / ν), χ_L_20 * 20^(- γ / ν), onlyMarks = true, mark = "o"),
        Plots.Linear(ts * 30^(1 / ν), χ_L_30 * 30^(- γ / ν), onlyMarks = true, mark = "o"),
        Plots.Linear(ts * 40^(1 / ν), χ_L_40 * 40^(- γ / ν), onlyMarks = true, mark = "o"),
        Plots.Linear(ts * 50^(1 / ν), χ_L_50 * 50^(- γ / ν), onlyMarks = true, mark = "o"),
    ],
    xlabel = L"t L^{1 / \nu}",
    ylabel = L"\chi L^{- \gamma / \nu}",
    style = "axis background/.style={fill=white}"
)

working_path = "D:\\Projects\\Modern Physics Experiments\\HTML5\\ising-data-collapsing\\"
name = "ising-data-collapsing-run-2.pdf"
save(working_path * name, p)