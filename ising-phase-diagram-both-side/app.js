const T_min = 0.001;
const T_max = 4.001;
const h_min = -2.0;
const h_max = 2.0;
const T_sampling_num = 100;
const h_sampling_num = 100;

const side_length = 10;
const heating_time = 1000;
const sampling_time = 3000;

const delta_T = (T_max - T_min) / T_sampling_num;
const delta_h = (h_max - h_min) / h_sampling_num;

const T_range = new Array(T_sampling_num).fill().map((_, i) => T_min + delta_T * i);
const h_range = new Array(h_sampling_num).fill().map((_, i) => h_min + delta_h * i);

const {site_list, inverse_list, neighbor_list} = square_lattice_2D_pbc(side_length);

let T_idx = 0, h_idx = 0;
let pause = false;

let canvas = document.getElementById("myCanvas"),
    context = canvas.getContext("2d");

function resetDisplay() {
    T_idx = 0;
    h_idx = 0;
    // Canvas background color being completely black 
    context.save();
    context.fillStyle = "rgba(0,0,0,1.0)";
    context.fillRect(0, 0, canvas.width, canvas.height);
    context.restore();
}

const single_point_width = canvas.width / T_sampling_num;
const single_point_height = canvas.height / h_sampling_num;

resetDisplay();

const r_value_0 = 0;
const g_value_0 = 0;
const b_value_0 = 255;
const r_value_1 = 0;
const g_value_1 = 255;
const b_value_1 = 255;
const r_value_2 = 255;
const g_value_2 = 0;
const b_value_2 = 255;

function magnetization_to_color(magnetization) {
    if (magnetization >= 0) {
        const r_value_current = r_value_0 + magnetization * (r_value_1 - r_value_0);
        const g_value_current = g_value_0 + magnetization * (g_value_1 - g_value_0);
        const b_value_current = b_value_0 + magnetization * (b_value_1 - b_value_0);
        return { r_value_current, g_value_current, b_value_current };
    } else {
        const r_value_current = r_value_0 - magnetization * (r_value_2 - r_value_0);
        const g_value_current = g_value_0 - magnetization * (g_value_2 - g_value_0);
        const b_value_current = b_value_0 - magnetization * (b_value_2 - b_value_0);
        return { r_value_current, g_value_current, b_value_current };
    }
}

// Animation: plotting the phase diagram
function plotPhaseDiagramPoint() {
    if (! pause) {
        const T = T_min + delta_T * T_idx;
        const beta = 1 / T;
        const h = h_min + delta_h * h_idx;

        // Heating up
        ising_config = new MPIsing(- beta, beta * h, site_list, inverse_list, neighbor_list);
        for (let heat_count = 0; heat_count < heating_time; heat_count++) {
            ising_config.update();
        }

        // Record the magnetization data
        magnetization_results = new Array(sampling_time).fill();
        for (let sample_count = 0; sample_count < sampling_time; sample_count++) {
            ising_config.update();
            magnetization_results[sample_count] = ising_config.magnetization();
        }
        const magnetization_avg = sum(magnetization_results) / sampling_time;

        // Paint the corresponding box in the canvas according to the MC simulation result
        context.save();
        // m=1 corresponds to blue, and m=0 corresponds to yellow
        const { r_value_current, g_value_current, b_value_current } = magnetization_to_color(magnetization_avg);
        context.fillStyle = `rgba(${r_value_current},${g_value_current},${b_value_current},1.0)`;
        context.fillRect(
            T_idx * single_point_width, 
            canvas.height - (h_idx + 1) * single_point_height,
            single_point_width, single_point_height
        );
        context.restore();

        // console.log(`Complete computation of T = ${T}, h = ${h}, average magnetization = ${magnetization_avg}`);

        h_idx ++;
        if (h_idx == h_range.length) {
            h_idx = 0;
            T_idx ++;
        }
        if (T_idx == T_range.length) {
            pause = true;
        }
        requestAnimationFrame(plotPhaseDiagramPoint);
    }
}

const start_phase_diagram = document.getElementById("start_phase_diagram");

start_phase_diagram.addEventListener("click", (ev) => {
    // cancelAnimationFrame(plotPhaseDiagramPoint);
    resetDisplay();
    plotPhaseDiagramPoint();
});