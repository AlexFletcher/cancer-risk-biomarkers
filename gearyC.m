% This function computes Geary's C

function gearyC=gearyC (cell_lattice, weights_matrix, total_weight_sum)
    N_size = length(cell_lattice(:,1));
    lattice = cellfun(@sum, cell_lattice); % Computes the total mutational burden of each of the cells on the lattice
    
    p =  reshape(lattice, 1,[])';
    p1 = p;

    q=sum(sum(((bsxfun(@minus,p',p)) .^2 ).* weights_matrix ));
    
    temp=(p1-mean(p1));
    minus_avg = sum(temp.*temp);
    gearyC = (N_size*N_size-1) * q / (minus_avg * 2 * total_weight_sum);

end