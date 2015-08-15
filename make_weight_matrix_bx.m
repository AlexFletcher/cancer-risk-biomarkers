
N_bx_radius =10; % set this variable to the radius of the biopsy
fName = 'weight_matrix_bx_small.mat' % set the filename of the weight matrix file here
numpoints = 2*N_bx_radius+ 1;

 xmat = bsxfun(@plus,-N_bx_radius:N_bx_radius,zeros(numpoints,numpoints))';
 ymat = xmat';
 dist_mat = sqrt(xmat.^2+ymat.^2);
 dist_mat = (dist_mat <= N_bx_radius);
 dist_mat =reshape(dist_mat,numpoints.*numpoints,1);
 xmat = reshape(xmat,numpoints.*numpoints,1);
 ymat = reshape(ymat,numpoints.*numpoints,1);
 
 xmat = bsxfun(@minus,xmat,xmat');
 ymat = bsxfun(@minus,ymat,ymat');
 
 weights_matrix_bx = 1 ./ (sqrt(xmat.^2 + ymat.^2) + 1);
 
 weights_matrix_bx = bsxfun(@times, weights_matrix_bx,dist_mat);
 weights_matrix_bx = bsxfun(@times, weights_matrix_bx,dist_mat');
 total_weight_sum_bx = sum(sum(weights_matrix_bx));
 

save(fName,'weights_matrix_bx','total_weight_sum_bx')