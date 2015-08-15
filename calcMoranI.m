% This function computes the value of Moran's I, for the given lattice and associated weight matrix

function moranI = calcMoranI (lattice_cells,weights_matrix,total_weight_sum)
    N_size = length(lattice_cells(:,1));

    lattice = zeros(N_size, N_size);
    for i=1:N_size
        for j=1:N_size
            all_muts= lattice_cells{i,j};
            lattice (i,j) = sum(all_muts); % here we compare cells based on the sum of their total number of mutations.
        end
    end

    avg_val  = mean(mean(lattice));

    lattice_col = reshape(lattice, 1,[])';
    lattice_col = lattice_col - avg_val;

    %load('weight_matrix.mat')
    moranI = lattice_col' * weights_matrix * lattice_col;
    moranI = moranI / (total_weight_sum * (lattice_col' * lattice_col));
    moranI=  moranI .* N_size .* N_size;

end
