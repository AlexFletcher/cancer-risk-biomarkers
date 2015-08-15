N_size =100; % size of lattice

weights_matrix = zeros(N_size*N_size,N_size*N_size);
nxt_indx_x = 1;
nxt_indx_y = 1;

total_weight_sum_bx = 0;
for i = 1:N_size
    i
    for j =1:N_size
        for k =1:N_size
            for l = 1:N_size
                p1 = [i j];
                p2 = [k l];
                temp = 1/(norm(p1-p2)+1);
                weights_matrix (nxt_indx_x, nxt_indx_y)=temp;
                total_weight_sum = total_weight_sum_bx + temp;
                nxt_indx_y = nxt_indx_y+1;
            end
        end
        nxt_indx_y = 1;
        nxt_indx_x = nxt_indx_x+1;
    end
    
end
save('weight_matrix.mat','weights_matrix','total_weight_sum')