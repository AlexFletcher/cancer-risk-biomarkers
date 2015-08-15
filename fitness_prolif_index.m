% This function computes the value of the proliferation-based indices

function ind_val = fitness_prolif_index(cell_matrix, prolif_locations, which_locs)
    N_size = length(cell_matrix(:,1));
    ind_val = 0;
    wt_sum = 0;
    x_mat = zeros(N_size,N_size);
    for i=1:N_size
        x_mat(i,:) = ones(1,N_size).*i;
    end
    y_mat = x_mat';

    for k=1:length(prolif_locations(:,1))
        if (~isequal(prolif_locations(k,:),[0 0]))
            x = x_mat - prolif_locations(k,1);
            y = y_mat - prolif_locations(k,2);

            dist = sqrt((x.*x) + (y.*y));

            dist = 1./(dist+1);
            dist = dist.*which_locs;
        
            temp_cell_mat = cell_matrix;
            temp_cell_mat(prolif_locations(k,1), prolif_locations(k,2))=0;
            ind_val = ind_val+ sum(sum(dist.*temp_cell_mat));
            wt_sum = wt_sum + sum(sum(dist))-1;
        end
    end
    ind_val = ind_val /wt_sum;
end
