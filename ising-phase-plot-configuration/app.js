const real_Tc = 2.2691853142;

let T = 4
let h = 0.0
let beta = 1 / T;

const side_length = 60;
const heating_time = 1000;
const sampling_time = 3000;

const {site_list, inverse_list, neighbor_list} = square_lattice_2D_pbc(side_length);

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

resetDisplay();

const single_point_width = canvas.width / side_length;
const single_point_height = canvas.height / side_length;

const r_up = 255;
const g_up = 0;
const b_up = 0;
const r_down = 0;
const g_down = 0;
const b_down = 255;

const ising_config = new MPIsing(- beta, beta * h, site_list, inverse_list, neighbor_list);

const waiting_time = 200;

setInterval(() => {
    if (! pause) {

        for (let site_idx = 0; site_idx < site_list.length; site_idx++) {
            const spin = ising_config.spin(site_idx);
            const [x_coord, y_coord] = site_list[site_idx];

            // Paint the corresponding box in the canvas according to the current field configuration
            context.save(); 

            let r_value_current = spin == 1? r_up : r_down;
            let g_value_current = spin == 1? g_up : g_down;
            let b_value_current = spin == 1? b_up : b_down;
            context.fillStyle = `rgba(${r_value_current},${g_value_current},${b_value_current},1.0)`;

            context.fillRect(
                x_coord * single_point_width, 
                canvas.height - (y_coord + 1) * single_point_height,
                single_point_width, single_point_height
            );
            context.restore();
        }
        
        ising_config.update();
        requestAnimationFrame(plotConfig);
    }
}, waiting_time);

const start_show_config = document.getElementById("start_show_config");

start_show_config.addEventListener("click", (ev) => {
    resetDisplay();
    plotConfig();
});

// The actual animation gives me a sense that clusters are moving leftward and downward.
// This is not necessarily illusion, because it faithfully reflects how the configuration is updated.