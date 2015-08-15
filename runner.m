% These are all of the cases we run for the data in the paper.

load('weight_matrix.mat')
load('weight_matrix_bx.mat')

%summary_stats('run_data_Nm10_s02_mu01.mat', 0.2,0.2,0.1, 10, 20, weights_matrix, total_weight_sum, weights_matrix_bx, total_weight_sum_bx)
%summary_stats('run_data_Nm10_s02_mu001.mat', 0.2,0.2,0.01, 10,20, weights_matrix, total_weight_sum, weights_matrix_bx, total_weight_sum_bx)
summary_stats('run_data_Nm10_s02_mu005.mat', 0.2,0.2,0.005, 10,20, weights_matrix, total_weight_sum, weights_matrix_bx, total_weight_sum_bx)

%summary_stats('run_data_Nm5_s02_mu01.mat', 0.2,0.2,0.1, 5,20, weights_matrix, total_weight_sum, weights_matrix_bx, total_weight_sum_bx)
% summary_stats('run_data_Nm3_s02_mu01.mat', 0.2,0.2,0.1, 3,20, weights_matrix, total_weight_sum, weights_matrix_bx, total_weight_sum_bx)
% summary_stats('run_data_Nm15_s02_mu01.mat', 0.2,0.2,0.1, 15, 20, weights_matrix, total_weight_sum, weights_matrix_bx, total_weight_sum_bx)
% 
% summary_stats('run_data_Nm10_s002_mu01.mat', 0.02,0.02,0.1, 10, 20, weights_matrix, total_weight_sum, weights_matrix_bx, total_weight_sum_bx)
% summary_stats('run_data_Nm10_s0002_mu01.mat', 0.002,0.002,0.1, 10, 20, weights_matrix, total_weight_sum, weights_matrix_bx, total_weight_sum_bx)
% summary_stats('run_data_Nm10_s0_mu01.mat', 0,0,0.1, 10, 20, weights_matrix, total_weight_sum, weights_matrix_bx, total_weight_sum_bx)
% 
% % for testing the effect of biopsy size
% 
% load('weight_matrix_bx_small.mat')
% summary_stats('run_data_Nm10_s02_mu01_bx_small.mat', 0.2,0.2,0.1, 10, 5,weights_matrix, total_weight_sum, weights_matrix_bx, total_weight_sum_bx)
% 
% load('weight_matrix_bx_large.mat')
% summary_stats('run_data_Nm10_s02_mu01_bx_large.mat', 0.2,0.2,0.1, 10, 40,weights_matrix, total_weight_sum, weights_matrix_bx, total_weight_sum_bx)