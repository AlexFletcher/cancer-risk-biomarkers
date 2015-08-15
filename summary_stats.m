% This code runs, collates, and saves the data from the lattice simulations N_runs times.

function summary_stats(fName, s_pos, s_del, mu, N_muts_necessary, bx_size, weights_matrix, total_weight_sum, weights_matrix_bx, total_weight_sum_bx)

N_runs = 200; % Change this for more or less runs

mut_probs_tested = mu;

tEnd=  600; % Stop saving data at time = 600.
time_step = 1;

times_vec = time_step:time_step:tEnd ;
N_times = length(times_vec);

shannon_indices_runs = cell(N_runs,length(mut_probs_tested));
simpson_indices_runs =  cell(N_runs,length(mut_probs_tested));
moranI_time_series_runs =  cell(N_runs,length(mut_probs_tested));
gearyC_time_series_runs = cell(N_runs,length(mut_probs_tested));

fpi0_vals = cell(N_runs,length(mut_probs_tested));
fpi_vals = cell(N_runs,length(mut_probs_tested));
fpi2_vals = cell(N_runs,length(mut_probs_tested));
fpi3_vals = cell(N_runs,length(mut_probs_tested));

mpi_pos_vals = cell(N_runs,length(mut_probs_tested));
mpi_tot_vals = cell(N_runs,length(mut_probs_tested));

lattice_erpos_runs = cell(N_runs,length(mut_probs_tested));
lattice_ki67_runs = cell(N_runs,length(mut_probs_tested));

time_of_cancer_occurrence_measure1_runs =  cell(N_runs,length(mut_probs_tested));
time_of_cancer_occurrence_measure2_runs = cell(N_runs,length(mut_probs_tested));

biopsy_types_of_mutants_runs =  cell(N_runs,length(mut_probs_tested));
biopsy_shannon_indices_runs =  cell(N_runs,length(mut_probs_tested));
biopsy_simpson_indices_runs =  cell(N_runs,length(mut_probs_tested));
biopsy_moranI_runs =  cell(N_runs,length(mut_probs_tested));
biopsy_gearyC_runs =  cell(N_runs,length(mut_probs_tested));

fpi0_bx_vals = cell(N_runs,length(mut_probs_tested));
fpi_bx_vals = cell(N_runs,length(mut_probs_tested));
fpi2_bx_vals = cell(N_runs,length(mut_probs_tested));
fpi3_bx_vals = cell(N_runs,length(mut_probs_tested));
mpi_pos_bx_vals = cell(N_runs,length(mut_probs_tested));
mpi_tot_bx_vals = cell(N_runs,length(mut_probs_tested));

biopsy_erpos_runs = cell(N_runs,length(mut_probs_tested));
biopsy_ki67_runs = cell(N_runs,length(mut_probs_tested));

scraping_types_of_mutants_runs =  cell(N_runs,length(mut_probs_tested));
scraping_shannon_indices_runs=  cell(N_runs,length(mut_probs_tested));
scraping_simpson_indices_runs =  cell(N_runs,length(mut_probs_tested));
scraping_erpos_runs = cell(N_runs,length(mut_probs_tested));
scraping_ki67_runs = cell(N_runs,length(mut_probs_tested));

types_of_mutants_time_series_runs =  cell(N_runs,length(mut_probs_tested));

for j=1:length(mut_probs_tested)
    
    parfor i=1:N_runs % Parallel loop
        
        i % Just to keep track

        [shannon_indices, simpson_indices, moranI_time_series,  gearyC_time_series, fpi0, fpi, fpi2, fpi3,mpi_pos,mpi_tot, lattice_erpos_time_series, lattice_ki67_time_series,...
            biopsy_shannon_indices, biopsy_simpson_indices, biopsy_moranI_series, biopsy_gearyC_series, fpi0_bx, fpi_bx, fpi2_bx, fpi3_bx,mpi_pos_bx,mpi_tot_bx,  biopsy_erpos_time_series, biopsy_ki67_time_series,...
            scraping_shannon_indices, scraping_simpson_indices, scraping_erpos_time_series, scraping_ki67_time_series,...
            time_of_cancer_occurrence_measure1, time_of_cancer_occurrence_measure2, types_of_mutants_time_series, biopsy_types_of_mutants, scraping_types_of_mutants]=...
            new_model2d(N_muts_necessary,s_pos, s_del,mut_probs_tested(j), bx_size, times_vec,weights_matrix,total_weight_sum,weights_matrix_bx,total_weight_sum_bx);
        
        shannon_indices_runs{i,j} = shannon_indices;
        simpson_indices_runs{i,j} =  simpson_indices;
        moranI_time_series_runs{i,j} = moranI_time_series;
        gearyC_time_series_runs{i,j} = gearyC_time_series;
        
        fpi0_vals{i,j} = fpi0;
        fpi_vals{i,j} = fpi;
        fpi2_vals{i,j} = fpi2;
        fpi3_vals{i,j} = fpi3;
        
        mpi_pos_vals{i,j} = mpi_pos;
        mpi_tot_vals{i,j} = mpi_tot;

        lattice_erpos_runs{i,j} = lattice_erpos_time_series;
        lattice_ki67_runs{i,j} = lattice_ki67_time_series;
        
        time_of_cancer_occurrence_measure1_runs{i,j} =  time_of_cancer_occurrence_measure1;
        time_of_cancer_occurrence_measure2_runs{i,j} = time_of_cancer_occurrence_measure2;
        
        biopsy_types_of_mutants_runs{i,j} =  biopsy_types_of_mutants;
        biopsy_shannon_indices_runs{i,j} =  biopsy_shannon_indices;
        biopsy_simpson_indices_runs{i,j} =  biopsy_simpson_indices;
        biopsy_moranI_runs{i,j} = biopsy_moranI_series;
        biopsy_gearyC_runs{i,j} = biopsy_gearyC_series;
        
        fpi0_bx_vals{i,j} = fpi0_bx;
        fpi_bx_vals{i,j} = fpi_bx;
        fpi2_bx_vals{i,j} = fpi2_bx;
        fpi3_bx_vals{i,j} = fpi3_bx;
        
        mpi_pos_bx_vals{i,j} = mpi_pos_bx;
        mpi_tot_bx_vals{i,j} = mpi_tot_bx;
        
        biopsy_erpos_runs{i,j} = biopsy_erpos_time_series;
        biopsy_ki67_runs{i,j} = biopsy_ki67_time_series;
        
        scraping_types_of_mutants_runs{i,j} =  scraping_types_of_mutants;
        scraping_shannon_indices_runs{i,j}=  scraping_shannon_indices;
        scraping_simpson_indices_runs{i,j} = scraping_simpson_indices;
        scraping_erpos_runs{i,j} = scraping_erpos_time_series;
        scraping_ki67_runs{i,j} = scraping_ki67_time_series;
        
        types_of_mutants_time_series_runs{i,j} = types_of_mutants_time_series;
        
    end
end
save(fName, 'shannon_indices_runs','simpson_indices_runs','moranI_time_series_runs','gearyC_time_series_runs','fpi0_vals','fpi_vals','fpi2_vals','fpi3_vals','mpi_pos_vals','mpi_tot_vals',...
    'lattice_erpos_runs','lattice_ki67_runs','time_of_cancer_occurrence_measure1_runs','time_of_cancer_occurrence_measure2_runs',...
    'biopsy_shannon_indices_runs','biopsy_simpson_indices_runs','biopsy_moranI_runs','biopsy_gearyC_runs','fpi0_bx_vals','fpi_bx_vals','fpi2_bx_vals','fpi3_bx_vals','mpi_pos_bx_vals','mpi_tot_bx_vals',...
    'biopsy_erpos_runs','biopsy_ki67_runs','scraping_shannon_indices_runs','scraping_simpson_indices_runs','scraping_erpos_runs','scraping_ki67_runs','types_of_mutants_time_series_runs','biopsy_types_of_mutants_runs','scraping_types_of_mutants_runs');

end
