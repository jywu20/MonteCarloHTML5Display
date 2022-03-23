/**
 * Simulating vanilla Ising model
 */

/**
 * Metropolis method for Ising model
 */
class MPIsing {
    constructor(J, h, site_list, inverse_list, neighbor_list) {
        this.J = J;
        this.h = h;
        if (site_list.length != inverse_list.length || 
            site_list.length != neighbor_list.length ||
            inverse_list.length != neighbor_list.length) {
            throw SyntaxError("The length of site_list, inverse_list and neighbor_list must be all the same.");
        }
        this.site_list = site_list;
        this.inverse_list = inverse_list;
        this.neighbor_list = neighbor_list;

        // All spin up by default
        this.configuration = new Array(site_list.length).fill(1);
    }

    update() {
        const J = this.J;
        const h = this.h;
        for (let current_site = 0; current_site < this.site_list.length; current_site++) {
            const sigma_i = this.configuration[current_site];

            let delta_F_J = 0.0;
            for (const neighbor_site of this.neighbor_list[current_site]) {
                const sigma_j = this.configuration[neighbor_site];
                delta_F_J += 2 * J * sigma_i * sigma_j;
            }

            const delta_F_h = 2 * h * sigma_i;
            const accept_rate = Math.exp(delta_F_J + delta_F_h);

            if (Math.random() < accept_rate) {
                this.configuration[current_site] *= -1;
            }
        }
    }

    free_energy() {
        let L_J = 0.0;
        let L_h = 0.0;
        for (let current_site = 0; current_site < this.site_list.length; current_site++) {
            const sigma_i = this.configuration[current_site];

            L_h += this.h * sigma_i;

            for (const neighbor_site of this.neighbor_list[current_site]) {
                const sigma_j = this.configuration[neighbor_site];
                L_J += this.J * sigma_i * sigma_j;
            }
        }
        L_J /= 2;

        return L_J + L_h;
    }

    magnetization() {
        return sum(this.configuration);
    }
}
