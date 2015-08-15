function [shannon_indices, simpson_indices, moranI_time_series,  gearyC_time_series, fpi0, fpi, fpi2, fpi3, mpi_pos, mpi_tot, lattice_erpos_time_series, lattice_ki67_time_series,...
    biopsy_shannon_indices, biopsy_simpson_indices, biopsy_moranI_series, biopsy_gearyC_series, fpi0_bx, fpi_bx, fpi2_bx, fpi3_bx,mpi_pos_bx, mpi_tot_bx, biopsy_erpos_time_series, biopsy_ki67_time_series,...
    scraping_shannon_indices, scraping_simpson_indices, scraping_erpos_time_series, scraping_ki67_time_series,...
    time_of_cancer_occurrence_measure1, time_of_cancer_occurrence_measure2, types_of_mutants_time_series, biopsy_type_of_mutants, scraping_cells_prop]=...
    new_model2d_birth_death(N_pos_mutations_for_cancer, s_pos, s_del, mut_prob, bx_size, time_series_to_save, weights_matrix, total_weight_sum, weights_matrix_bx, total_weight_sum_bx)

t_cutoff_prolif = 0.01; % this is the time window for which cells undergoing division stain positive for Ki-67
N_cell_cutoff_for_cancer = 1;
N_pos_mut_cutoff_for_er_pos = 1; % This is the number of mutations, above which, the cell is considered "traditional biomarker" positive, referred in the code as ER positive
N_percentage_cutoff_for_cancer = 0.05; % 5% cancer cells = cancer
time_of_cancer_occurrence_measure1 =0; % this will be the measure such that we have a single cell with cancer

time_of_cancer_occurrence_measure2 =0; % this is the time at which we have a certain percentatge of cancer cells -- the primary endpoint in the paper

num_cancer_cells = 0; % this variable to keep track of the number of cancer cells

cur_iter = 0; % tracks iteration number of the loop

N_timepoints_saved = length(time_series_to_save);
next_time = 1;
shannon_indices = zeros(N_timepoints_saved,1);
simpson_indices = zeros(N_timepoints_saved,1);
lattice_ki67_time_series = zeros(N_timepoints_saved,1);
scraping_ki67_time_series =  zeros(N_timepoints_saved,1);
lattice_erpos_time_series = zeros(N_timepoints_saved,1);
scraping_erpos_time_series =  zeros(N_timepoints_saved,1);
types_of_mutants_time_series = cell(N_timepoints_saved,1);

moranI_time_series = zeros(N_timepoints_saved,1);
gearyC_time_series = zeros(N_timepoints_saved,1);
cell_props=cell(N_timepoints_saved,1);
N_cells_scraped = 1000;
scraping_cells = cell(N_cells_scraped,N_timepoints_saved);
scraping_cells_prop=cell(N_timepoints_saved,1);
scraping_shannon_indices = zeros(N_timepoints_saved,1);
scraping_simpson_indices = zeros(N_timepoints_saved,1);

N_size = 100;
cur_time = 0;

load('initial_locs.mat')

cell_lattice = cell(N_size, N_size);

for i=1:N_size
    for j=1:N_size
        cell_lattice{i,j} = [0 0 0];
    end
end

times_of_cellular_events = zeros(10000000,3);
%%%%%%%%%%%%%%%% IMPORTANT Biopsy variables %%%%%%%%%%%%%%%%%%%
% For multiple biopsies at a given time, set this here by uncommenting the next two lines.
%N_extra_bx_at_t_50 = 10;
%biopsy_time_points =[time_series_to_save(1:50) (50.*(ones(1,N_extra_bx_at_t_50-1))) time_series_to_save(51:end)]; %[2 2 3 3];

% For the standard case of just one biopsy at all timepoints, then uncomment the following line:
biopsy_time_points = time_series_to_save;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

biopsy_radius = bx_size;
num_biopsies = length(biopsy_time_points);
next_biopsy_point = 1;
biopsy_region = cell(2*biopsy_radius+1,2*biopsy_radius+1,num_biopsies);
biopsy_type_of_mutants = cell(num_biopsies,1);
biopsy_shannon_indices = zeros(num_biopsies,1);
biopsy_simpson_indices = zeros(num_biopsies,1);
biopsy_moranI_series = zeros(num_biopsies,1);
biopsy_gearyC_series = zeros(num_biopsies,1);
fpi0 = zeros(N_timepoints_saved,1);
fpi = zeros(N_timepoints_saved,1);
fpi2 = zeros(N_timepoints_saved,1);
fpi3 = zeros(N_timepoints_saved,1);

mpi_pos = zeros(N_timepoints_saved,1);
mpi_tot = zeros(N_timepoints_saved,1);

fpi0_bx = zeros(num_biopsies,1);
fpi_bx = zeros(num_biopsies,1);
fpi2_bx = zeros(num_biopsies,1);
fpi3_bx = zeros(num_biopsies,1);

mpi_pos_bx = zeros(num_biopsies,1);
mpi_tot_bx = zeros(num_biopsies,1);

biopsy_ki67_time_series =  zeros(num_biopsies,1);
biopsy_erpos_time_series =  zeros(num_biopsies,1);

types_of_mutants = [];
types_of_mutants(1,:) = [0 0 0 (N_size*N_size)]; % keep this in the format d n p #of cells; d, n, p, are the numbers of each type of mutation
cell_locations = cell(1,1,1); % so far we only have type 0 0 0 (all unmutated)
cell_locations{1,1,1} = initial_locs; % locations of the cells are initialized
cell_locations{1,1,2} = [];
propen = [];
propen(1) = types_of_mutants(1,4) * ((1 - s_del)^(types_of_mutants(1,1))) * ((1 + s_pos)^(types_of_mutants(1,3)));

while (1) % we stop when the endpoint of cancer is reached
    % First decide if we have to save this timepoint or take a biopsy
    
    if (next_time <= N_timepoints_saved)
        if(cur_time >= time_series_to_save(next_time))
            
            cell_matrix_output = zeros(N_size, N_size);
            cell_matrix_output_2 = zeros(N_size, N_size);
            cell_matrix_output_3 = zeros(N_size, N_size);
            cell_matrix_output_0 = zeros(N_size, N_size);

            mpi_pos_matrix_output = zeros(N_size, N_size);
            mpi_tot_matrix_output = zeros(N_size, N_size);

            which_locs = ones(N_size, N_size); % means that all cells of lattice are involved in FPI calculation

            for i=1:N_size
                for j=1:N_size
                    all_muts= cell_lattice{i,j};
                    cell_matrix_output_0 (i,j) = (all_muts(3)>=N_pos_mutations_for_cancer-1); %sum(all_muts);
                    cell_matrix_output (i,j) = (all_muts(3)>=N_pos_mutations_for_cancer-2); %sum(all_muts);
                    cell_matrix_output_2 (i,j) = (all_muts(3)>=N_pos_mutations_for_cancer-3); %sum(all_muts);
                    cell_matrix_output_3 (i,j) = (all_muts(3) + all_muts(2) >= fix(N_pos_mutations_for_cancer/2)); %sum(all_muts);
                    mpi_pos_matrix_output(i,j) = all_muts(3);
                    mpi_tot_matrix_output(i,j) = sum(all_muts);
                end
            end
            
            for i=cur_iter:-1:1
                if ((cur_time-times_of_cellular_events(i,1))>t_cutoff_prolif)
                    last_prolif = times_of_cellular_events((i+1):cur_iter,2:3);
                    break;
                end
            end
            
            lattice_ki67_time_series(next_time) = length(last_prolif(:,1))/(N_size*N_size);
            lattice_erpos_time_series(next_time) = sum(types_of_mutants(:,4).*(types_of_mutants(:,3)>N_pos_mut_cutoff_for_er_pos))/(N_size*N_size);
            fpi0(next_time) = fitness_prolif_index(cell_matrix_output_0, last_prolif,which_locs);
            fpi(next_time) = fitness_prolif_index(cell_matrix_output, last_prolif,which_locs);
            fpi2(next_time) = fitness_prolif_index(cell_matrix_output_2, last_prolif,which_locs);
            fpi3(next_time) = fitness_prolif_index(cell_matrix_output_3, last_prolif,which_locs);
            
            mpi_pos(next_time) = fitness_prolif_index(mpi_pos_matrix_output, last_prolif,which_locs);
            mpi_tot(next_time) = fitness_prolif_index(mpi_tot_matrix_output, last_prolif,which_locs);
            
            types_of_mutants_time_series{next_time} = types_of_mutants;
            
            %%cell_lattice_time_series (:,:,next_time) = cell_lattice; % if you wanted to save the cell lattice itself as a time series, you would uncomment this line

            moranI_time_series(next_time) = calcMoranI(cell_lattice,weights_matrix,total_weight_sum);
            gearyC_time_series(next_time) = gearyC(cell_lattice, weights_matrix,total_weight_sum);

            cur_cell_props = types_of_mutants(:,4)./(N_size*N_size);
            cell_props{next_time} = cur_cell_props;
            [temp1, temp2] = diversityMeasures(cur_cell_props);

            shannon_indices(next_time) = temp1;
            simpson_indices(next_time) = temp2;
            
            % do the same calculations but for scraped cells
            % first get the cells that were scraped
            
            scraped_indices=randperm(N_size*N_size,N_cells_scraped) - 1;
            scraped_index_x =mod(scraped_indices,N_size) + 1;
            scraped_index_y = fix(scraped_indices/N_size) + 1;
            types_of_mutants_scraping=[];
            
            for r=1:N_cells_scraped
                scraping_cells{r, next_time} = cell_lattice{int16(scraped_index_x(r)), int16(scraped_index_y(r))};

                if(r==1)
                    types_of_mutants_scraping(1,:) = [scraping_cells{r,next_time} 1];
                else
                    foundflag=false;
                    for s=1:length(types_of_mutants_scraping(:,1))
                        if (isequal(types_of_mutants_scraping(s,1:3),scraping_cells{r,next_time}))
                            foundflag=true;
                            types_of_mutants_scraping(s,4) = types_of_mutants_scraping(s,4)+1;
                            break;
                        end
                        
                    end
                    if(~foundflag)
                        types_of_mutants_scraping = [types_of_mutants_scraping; scraping_cells{r,next_time} 1];
                    end
                end
                
            end
            
            prolif_count_scraped = sum(ismember(last_prolif, [int16(scraped_index_x') int16(scraped_index_y')] ,'rows'));
            scraping_ki67_time_series(next_time) = prolif_count_scraped/N_cells_scraped;
            
            scraping_cells_prop{next_time} = types_of_mutants_scraping;% types_of_mutants_scraping(:,4)./sum(types_of_mutants_scraping(:,4));
            
            [temp1, temp2] = diversityMeasures( types_of_mutants_scraping(:,4)./sum(types_of_mutants_scraping(:,4)));
            
            scraping_erpos_time_series(next_time) = sum(types_of_mutants_scraping(:,4).*(types_of_mutants_scraping(:,3)>N_pos_mut_cutoff_for_er_pos))/N_cells_scraped;
            
            scraping_shannon_indices(next_time) = temp1;
            scraping_simpson_indices(next_time) = temp2;
            
            next_time = next_time + 1;
            
        end
    end
    
    
    if(next_biopsy_point <= num_biopsies)
        if(cur_time >= biopsy_time_points(next_biopsy_point))
            next_pt_diff = false;
            while(~next_pt_diff)
                biopsy_centre_x = randi(N_size- 2*biopsy_radius ) + (biopsy_radius);
                biopsy_centre_y = randi(N_size- 2*biopsy_radius ) + (biopsy_radius);
                types_of_mutants_biopsy=[];
                temp0=cell_lattice{(biopsy_centre_x ), (biopsy_centre_y)};
                types_of_mutants_biopsy(1,:) = [temp0 1];
                for i=0:biopsy_radius
                    for j=0:biopsy_radius
                        if (i^2 + j^2 <= biopsy_radius^2)
                            temp1=cell_lattice{(biopsy_centre_x + i), (biopsy_centre_y+j)};
                            biopsy_region{biopsy_radius + 1 + i,biopsy_radius + 1 + j,next_biopsy_point} = temp1;
                            temp2=cell_lattice{(biopsy_centre_x + i), (biopsy_centre_y-j)};
                            biopsy_region{biopsy_radius + 1 + i,biopsy_radius + 1 - j,next_biopsy_point} = temp2;
                            temp3=cell_lattice{(biopsy_centre_x - i), (biopsy_centre_y+j)};
                            biopsy_region{biopsy_radius + 1 - i,biopsy_radius + 1 + j,next_biopsy_point} = temp3;
                            temp4 = cell_lattice{(biopsy_centre_x - i), (biopsy_centre_y-j)};
                            biopsy_region{biopsy_radius + 1 - i,biopsy_radius + 1 - j,next_biopsy_point} = temp4;
                            
                            
                            if (i > 0 && j > 0)
                                foundflag= false;
                                
                                for r=1:length(types_of_mutants_biopsy(:,1))
                                    
                                    if(isequal(types_of_mutants_biopsy(r,1:3), temp1 ))
                                        foundflag= true;
                                        types_of_mutants_biopsy(r,4) = types_of_mutants_biopsy(4) + 1;
                                    end
                                end
                                if(~foundflag)
                                    types_of_mutants_biopsy = [types_of_mutants_biopsy; temp1 1];
                                end
                                foundflag= false;
                                
                                for r=1:length(types_of_mutants_biopsy(:,1))
                                    
                                    if(isequal(types_of_mutants_biopsy(r,1:3), temp2 ))
                                        foundflag= true;
                                        types_of_mutants_biopsy(r,4) = types_of_mutants_biopsy(4) + 1;
                                    end
                                end
                                if(~foundflag)
                                    types_of_mutants_biopsy = [types_of_mutants_biopsy; temp2 1];
                                end
                                
                                foundflag= false;
                                
                                for r=1:length(types_of_mutants_biopsy(:,1))
                                    
                                    if(isequal(types_of_mutants_biopsy(r,1:3), temp3 ))
                                        foundflag= true;
                                        types_of_mutants_biopsy(r,4) = types_of_mutants_biopsy(4) + 1;
                                    end
                                end
                                if(~foundflag)
                                    types_of_mutants_biopsy = [types_of_mutants_biopsy; temp3 1];
                                end
                                
                                foundflag= false;
                                
                                for r=1:length(types_of_mutants_biopsy(:,1))
                                    
                                    if(isequal(types_of_mutants_biopsy(r,1:3), temp4 ))
                                        foundflag= true;
                                        types_of_mutants_biopsy(r,4) = types_of_mutants_biopsy(4) + 1;
                                    end
                                end
                                if(~foundflag)
                                    types_of_mutants_biopsy = [types_of_mutants_biopsy; temp4 1];
                                end
                                
                            end
                            
                            
                        end
                    end
                end
                
                biopsy_type_of_mutants{next_biopsy_point} = types_of_mutants_biopsy;
                [temp1,temp2] = diversityMeasures(types_of_mutants_biopsy(:,4)./sum(types_of_mutants_biopsy(:,4)));
                biopsy_shannon_indices(next_biopsy_point) = temp1;
                biopsy_simpson_indices(next_biopsy_point) = temp2;

                bx_region = cell_lattice((biopsy_centre_x - biopsy_radius):(biopsy_centre_x + biopsy_radius),(biopsy_centre_y - biopsy_radius):(biopsy_centre_y + biopsy_radius));
                biopsy_moranI_series(next_time) = calcMoranI(bx_region,weights_matrix_bx,total_weight_sum_bx);
                biopsy_gearyC_series(next_time) = gearyC(bx_region,weights_matrix_bx,total_weight_sum_bx);
                
                biopsy_erpos_time_series(next_biopsy_point) = sum(types_of_mutants_biopsy(:,4).*(types_of_mutants_biopsy(:,3)>N_pos_mut_cutoff_for_er_pos))/(pi*biopsy_radius*biopsy_radius);
                
                cell_matrix_output = zeros(N_size, N_size);
                cell_matrix_output_2 = zeros(N_size, N_size);
                cell_matrix_output_3 = zeros(N_size, N_size);
                cell_matrix_output_0 = zeros(N_size, N_size);

                mpi_pos_matrix_output = zeros(N_size, N_size);
                mpi_tot_matrix_output = zeros(N_size, N_size);
                which_locs = zeros(N_size,N_size);
                for r=-biopsy_radius:biopsy_radius
                    for s=-biopsy_radius:biopsy_radius
                        if (norm([r s]) <= biopsy_radius)
                        all_muts= cell_lattice{r+biopsy_centre_x,s+biopsy_centre_y};
                        cell_matrix_output_0 (r+biopsy_centre_x,s+biopsy_centre_y) = (all_muts(3)>=N_pos_mutations_for_cancer-1); %sum(all_muts);
                        
                        cell_matrix_output (r+biopsy_centre_x,s+biopsy_centre_y) = (all_muts(3)>=N_pos_mutations_for_cancer-2); %sum(all_muts);
                        cell_matrix_output_2 (r+biopsy_centre_x,s+biopsy_centre_y) = (all_muts(3)>=N_pos_mutations_for_cancer-3); %sum(all_muts);
                        cell_matrix_output_3 (r+biopsy_centre_x,s+biopsy_centre_y) = (all_muts(3) + all_muts(2) >=fix(N_pos_mutations_for_cancer/2)); %sum(all_muts);
                        mpi_pos_matrix_output(r+biopsy_centre_x,s+biopsy_centre_y) = all_muts(3);
                        mpi_tot_matrix_output(r+biopsy_centre_x,s+biopsy_centre_y) = sum(all_muts);
                        which_locs(r+biopsy_centre_x,s+biopsy_centre_y)= 1;
                        end
                    end
                end

                bx_prolif = [];
                
                for p=cur_iter:-1:1
                    if ((cur_time-times_of_cellular_events(p,1))>t_cutoff_prolif)
                        last_prolif = times_of_cellular_events((p+1):cur_iter,2:3);
                        break;
                    end
                end
                for r = 1:length(last_prolif(:,1))
                    if (norm(last_prolif(r,:) - [biopsy_centre_x biopsy_centre_y]) <= biopsy_radius)
                        bx_prolif(length(bx_prolif)+1,:) = last_prolif(r,:);
                    end
                end
                
                
                if(~isempty(bx_prolif))
                    biopsy_ki67_time_series(next_biopsy_point) = length(bx_prolif(:,1))/(biopsy_radius*biopsy_radius*pi);
                    
                    
                    fpi0_bx(next_biopsy_point) = fitness_prolif_index(cell_matrix_output_0, bx_prolif,which_locs);
                    
                    fpi_bx(next_biopsy_point) = fitness_prolif_index(cell_matrix_output, bx_prolif,which_locs);
                    fpi2_bx(next_biopsy_point) = fitness_prolif_index(cell_matrix_output_2, bx_prolif,which_locs);
                    
                    fpi3_bx(next_biopsy_point) = fitness_prolif_index(cell_matrix_output_3, bx_prolif,which_locs);
                    mpi_pos_bx(next_biopsy_point) = fitness_prolif_index(mpi_pos_matrix_output, bx_prolif,which_locs);
                    mpi_tot_bx(next_biopsy_point) = fitness_prolif_index(mpi_tot_matrix_output, bx_prolif,which_locs);
                else
                    biopsy_ki67_time_series(next_biopsy_point) = 0;
                    
                    fpi0_bx(next_biopsy_point) = 0;
                    
                    fpi_bx(next_biopsy_point) = 0;
                    fpi2_bx(next_biopsy_point) = 0;
                    
                    fpi3_bx(next_biopsy_point) = 0;
                    mpi_pos_bx(next_biopsy_point) = 0;
                    mpi_tot_bx(next_biopsy_point) = 0; 
                end
                
                
                if(next_biopsy_point  < num_biopsies)
                    if(biopsy_time_points(next_biopsy_point) ~= biopsy_time_points(next_biopsy_point+1))
                        next_pt_diff=true;
                    end
                else
                    next_pt_diff=true;
                end
                next_biopsy_point = next_biopsy_point+1;
                
            end
        end
    end

% The following is the Gillespie algorithm:
    total_prop = sum(propen);
    rand_1 = rand(1);
    
    if (rand_1 < 1)
        rand_1 = rand_1 + eps;
    end
    rand_2 = rand(1) * total_prop;
    dT = -log(rand_1) / total_prop;
    if (isrow(propen))
        propen = propen';
    end;
    cum_propen = [0;cumsum(propen)];
    [~,ind]  = histc(rand_2,cum_propen);
    if (ind > length(cum_propen)-1)
        ind = ind-1;
    end
    
    % so now ind tells you which cell type has been chosen to die
    % we have to look through lists of location of cells of that type to
    % choose what cell should die.
    % first get those lists:

    mutant_chosen = types_of_mutants(ind,:);
    list_of_cell_locs_to_die = cell_locations{mutant_chosen(1) + 1,mutant_chosen(2) + 1,mutant_chosen(3) + 1};

    row_vals =list_of_cell_locs_to_die(:,1);
    
    row_index = randi(length(row_vals)); %this is for the spatially uniform death case.
    
    %% the following code is what tells the program to choose a cell more weighted towards the bottom of the lattice (if you wanted spatially weighted death, where bottom is 2x as likely as top, with linear gradient)
    %         baseline = 1 / (N_size-1);
    %         slope = (N_size - 2) / (N_size-1);
    %
    %         scaled_row_vals = baseline + (slope.*(row_vals));
    %
    %         cum_rowVals = [0;cumsum(scaled_row_vals)];
    %         [~,row_index]  = histc(rand,cum_rowVals ./ sum(scaled_row_vals));
    %         if (row_index > length(cum_rowVals)-1)
    %             row_index = row_index-1;
    %         end
    
    cell_loc_x = list_of_cell_locs_to_die(row_index,1);
    cell_loc_y = list_of_cell_locs_to_die(row_index,2);
    
    %list_of_cell_locs_to_die(row_index,:) = [];
    %types_of_mutants(ind,4) = types_of_mutants(ind,4)-1; % cell is officially killed.
   % propen(ind) = types_of_mutants(ind,4) * ((1 - s_del)^(types_of_mutants(ind,1))) * ((1 + s_pos)^(types_of_mutants(ind,3)));
    %cell_locations{mutant_chosen(1) + 1,mutant_chosen(2) + 1,mutant_chosen(3) + 1} = list_of_cell_locs_to_die;

    % now lets choose one of its neighbours at random
    
    %% This is the part where we choose a neighbour:
    if (cell_loc_x == 1 &&  cell_loc_y ==1)
        temp = randi(2);
        if (temp ==1)
            cell_neighbour_x = cell_loc_x + 1;
            cell_neighbour_y = cell_loc_y;
        else
            
            cell_neighbour_x = cell_loc_x ;
            cell_neighbour_y = cell_loc_y + 1;
        end
        
    elseif(cell_loc_x == N_size && cell_loc_y ==1)
        
        temp = randi(2);
        if (temp ==1)
            cell_neighbour_x = cell_loc_x - 1;
            cell_neighbour_y = cell_loc_y;
        else
            
            cell_neighbour_x = cell_loc_x ;
            cell_neighbour_y = cell_loc_y + 1;
        end
        
    elseif(cell_loc_x == 1 && cell_loc_y ==N_size)
        temp = randi(2);
        if (temp ==1)
            cell_neighbour_x = cell_loc_x +1;
            cell_neighbour_y = cell_loc_y;
        else
            
            cell_neighbour_x = cell_loc_x ;
            cell_neighbour_y = cell_loc_y - 1;
        end
    elseif(cell_loc_x == N_size && cell_loc_y ==N_size)
        temp = randi(2);
        if (temp ==1)
            cell_neighbour_x = cell_loc_x - 1;
            cell_neighbour_y = cell_loc_y;
        else
            
            cell_neighbour_x = cell_loc_x ;
            cell_neighbour_y = cell_loc_y - 1;
        end
    elseif(cell_loc_x == N_size )
        temp = randi(3);
        if (temp ==1)
            cell_neighbour_x = cell_loc_x ;
            cell_neighbour_y = cell_loc_y + 1;
        elseif (temp==2)
            cell_neighbour_x = cell_loc_x ;
            cell_neighbour_y = cell_loc_y - 1;
        else
            cell_neighbour_x = cell_loc_x  -1 ;
            cell_neighbour_y = cell_loc_y;
        end
    elseif(cell_loc_y == N_size)
        temp = randi(3);
        if (temp ==1)
            cell_neighbour_x = cell_loc_x +1  ;
            cell_neighbour_y = cell_loc_y ;
        elseif (temp==2)
            cell_neighbour_x = cell_loc_x -1 ;
            cell_neighbour_y = cell_loc_y ;
        else
            cell_neighbour_x = cell_loc_x  ;
            cell_neighbour_y = cell_loc_y - 1;
        end
    elseif(cell_loc_x ==1)
        temp = randi(3);
        if (temp ==1)
            cell_neighbour_x = cell_loc_x ;
            cell_neighbour_y = cell_loc_y + 1;
        elseif (temp==2)
            cell_neighbour_x = cell_loc_x ;
            cell_neighbour_y = cell_loc_y - 1;
        else
            cell_neighbour_x = cell_loc_x  +1 ;
            cell_neighbour_y = cell_loc_y;
        end
    elseif(cell_loc_y==1)
        temp = randi(3);
        if (temp ==1)
            cell_neighbour_x = cell_loc_x +1 ;
            cell_neighbour_y = cell_loc_y;
        elseif (temp==2)
            cell_neighbour_x = cell_loc_x -1 ;
            cell_neighbour_y = cell_loc_y;
        else
            cell_neighbour_x = cell_loc_x ;
            cell_neighbour_y = cell_loc_y + 1;
        end
    else
        temp = randi(4);
        if (temp ==1)
            cell_neighbour_x = cell_loc_x +1 ;
            cell_neighbour_y = cell_loc_y;
        elseif (temp==2)
            cell_neighbour_x = cell_loc_x -1 ;
            cell_neighbour_y = cell_loc_y;
        elseif (temp==3)
            cell_neighbour_x = cell_loc_x ;
            cell_neighbour_y = cell_loc_y+1;
        else
            cell_neighbour_x = cell_loc_x ;
            cell_neighbour_y = cell_loc_y -1;
        end
    end

    does_mutation_happen = (rand < mut_prob);
    what_mutation_happens = randi(3); % randomly choose between D, N, P types of mutations
    
    mutation_occurring = does_mutation_happen * what_mutation_happens;
    
    old_cell_type = cell_lattice{cell_loc_x, cell_loc_y};
    neighbour_cell_type = cell_lattice{cell_neighbour_x, cell_neighbour_y};
    %have to kill off neighbour cell
    for k=1:length(types_of_mutants(:,1))
        if (isequal(types_of_mutants(k,1:3),neighbour_cell_type))
            
            types_of_mutants(k,4) = types_of_mutants(k,4)-1;
            propen(k) = types_of_mutants(k,4) * ((1 - s_del)^(types_of_mutants(k,1))) * ((1 + s_pos)^(types_of_mutants(k,3)));
            
            
            break;
        end
        
    end
    
    if (neighbour_cell_type(3) >= N_pos_mutations_for_cancer)
                num_cancer_cells = num_cancer_cells -1;

    end
    
    new_locs_list = cell_locations{neighbour_cell_type(1) + 1, neighbour_cell_type(2) + 1, neighbour_cell_type(3) + 1};
    for k=1:length(new_locs_list(:,1))
        if (isequal(new_locs_list(k,:),[cell_neighbour_x cell_neighbour_y]))
            new_locs_list(k,:)=[]; 
            break;
        end
        
    end
    cell_locations{neighbour_cell_type(1) + 1, neighbour_cell_type(2) + 1, neighbour_cell_type(3) + 1} = new_locs_list;

    if(mutation_occurring  >=1 )
        % same as neighbour cell type, but instead has a deleterious
        % mutation added
        new_cell_type = cell_lattice{cell_loc_x, cell_loc_y};
        new_cell_type(mutation_occurring) = new_cell_type(mutation_occurring)+ 1; % adds a  mutation
        cell_lattice{cell_neighbour_x,cell_neighbour_y} = new_cell_type;
        
        foundflag= false;
        for k=1:length(types_of_mutants(:,1))
            if (isequal(types_of_mutants(k,1:3),new_cell_type))
                foundflag=true;
                
                types_of_mutants(k,:) = [new_cell_type (types_of_mutants(k,4) +1)];
                propen(k) = types_of_mutants(k,4) * ((1 - s_del)^(types_of_mutants(k,1))) * ((1 + s_pos)^(types_of_mutants(k,3)));
                
                new_locs_list = cell_locations{new_cell_type(1) + 1, new_cell_type(2) + 1, new_cell_type(3) + 1};
                cell_locations{new_cell_type(1) + 1, new_cell_type(2) + 1, new_cell_type(3) + 1} = [new_locs_list; cell_neighbour_x cell_neighbour_y];
                break;
            end
            
        end
        if (~foundflag)
            types_of_mutants = [types_of_mutants; new_cell_type 1];
            temp_var = length(types_of_mutants(:,1));
            propen(temp_var) = types_of_mutants(temp_var,4) * ((1 - s_del)^(types_of_mutants(temp_var,1))) * ((1 + s_pos)^(types_of_mutants(temp_var,3)));
            
            cell_locations{new_cell_type(1) + 1, new_cell_type(2) + 1, new_cell_type(3) + 1} = [cell_neighbour_x cell_neighbour_y];
            
        end
        
        if (new_cell_type(3) >=N_pos_mutations_for_cancer)
            num_cancer_cells = num_cancer_cells +1;
            if (num_cancer_cells >= N_cell_cutoff_for_cancer && time_of_cancer_occurrence_measure1 ==0)
                time_of_cancer_occurrence_measure1 = cur_time;
            end
            if(num_cancer_cells/(N_size*N_size)  >= N_percentage_cutoff_for_cancer && time_of_cancer_occurrence_measure2 ==0)
                time_of_cancer_occurrence_measure2 = cur_time;
                break;
            end
            
        end
        
    else
                new_cell_type = cell_lattice{cell_loc_x, cell_loc_y};

        cell_lattice{cell_neighbour_x,cell_neighbour_y} = new_cell_type;
        
        new_locs_list = cell_locations{new_cell_type(1) + 1, new_cell_type(2) + 1, new_cell_type(3) + 1};

        cell_locations{new_cell_type(1) + 1, new_cell_type(2) + 1, new_cell_type(3) + 1} = [new_locs_list; cell_neighbour_x cell_neighbour_y];

        for k=1:length(types_of_mutants(:,1))
            if (isequal(types_of_mutants(k,1:3), new_cell_type))
                types_of_mutants(k,4) = types_of_mutants(k,4) +1;
                propen(k) = types_of_mutants(k,4) * ((1 - s_del)^(types_of_mutants(k,1))) * ((1 + s_pos)^(types_of_mutants(k,3)));
                break;
            end
        end
        
        if (new_cell_type(3) >=N_pos_mutations_for_cancer)
            num_cancer_cells = num_cancer_cells +1;
            if (num_cancer_cells >= N_cell_cutoff_for_cancer && time_of_cancer_occurrence_measure1 ==0)
                time_of_cancer_occurrence_measure1 = cur_time;
            end
            if(num_cancer_cells/(N_size*N_size)  >= N_percentage_cutoff_for_cancer && time_of_cancer_occurrence_measure2 ==0)
                time_of_cancer_occurrence_measure2 = cur_time;
                break;
            end
            
        end
        
    end

    cur_time = cur_time + dT;
    cur_iter = cur_iter + 1;
    times_of_cellular_events(cur_iter,:) = [cur_time cell_loc_x cell_loc_y];
    
end

% The following are for outputting visualizations of the lattice

% cell_matrix_output = zeros(biopsy_radius, biopsy_radius);
% colormap(jet)
% 
% for k=1:num_biopsies
%     for i=1:2*biopsy_radius+1
% 
%         for j=1:2*biopsy_radius+1
% 
%             all_muts= biopsy_region{i,j,k}; 
%             if (~isempty(all_muts))
%             cell_matrix_output (i,j) = all_muts(3);%sum(all_muts);
%             end
%         end
%     end
% 
%     figure();
%     caxis([0, 10])
%     h = imagesc(cell_matrix_output);
%     colormap(jet)
% 
% end



%  cell_matrix_output = zeros(N_size, N_size);
%  for i=1:N_size
%     for j=1:N_size
%          all_muts= cell_lattice{i,j};
%          cell_matrix_output (i,j) = all_muts(3); %sum(all_muts);
% 
%      end
%  end
% 
% % fpi = fitness_prolif_index(cell_matrix_output, last_prolif);
% 
% types_of_mutants
% figure();
% h = imagesc(cell_matrix_output);
% caxis([0,10])
% colormap(jet)

end