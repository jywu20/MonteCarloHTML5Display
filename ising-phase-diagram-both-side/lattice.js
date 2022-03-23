function back_into_range(idx, range) {
    return (idx + range) % range;
}

function square_lattice_2D_pbc(side_length) {
    const site_num = side_length**2;
    
    const site_list = new Array(site_num).fill();
    const inverse_list = new Array(site_num).fill().map(_ => new Array(site_num).fill());
    const neighbor_list = new Array(site_num).fill().map(_ => new Array(4).fill());

    for (let site = 0; site < site_num; site++) {
        const y = site % side_length;
        const x = (site - y) / side_length;

        site_list[site] = [x, y];
        inverse_list[x][y] = site;
    }

    for (let site = 0; site < site_num; site++) {
        let [x, y] = site_list[site];
        neighbor_list[site][0] = inverse_list[back_into_range(x + 1, side_length)][y];
        neighbor_list[site][1] = inverse_list[back_into_range(x - 1, side_length)][y];
        neighbor_list[site][2] = inverse_list[x][back_into_range(y + 1, side_length)];
        neighbor_list[site][3] = inverse_list[x][back_into_range(y - 1, side_length)];
    }

    return {site_list, inverse_list, neighbor_list};
}