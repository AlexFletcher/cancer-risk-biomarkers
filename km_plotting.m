load('../Data/mat/run_data_Nm15_s02_mu01.mat') % load the appropriate data file here.

N_m = 10; %% Only have to set if plotting -- this is the variable for the number of mutations necessary to reach cancer

N_runs = length(shannon_indices_runs)

% The following is just variable/matrix formatting:

shannon_indices_runs = cell2mat(shannon_indices_runs')';
simpson_indices_runs = cell2mat(simpson_indices_runs')';
moranI_runs = cell2mat(moranI_time_series_runs')';
gearyC_runs = cell2mat(gearyC_time_series_runs')';

scraping_shannon_indices_runs = cell2mat(scraping_shannon_indices_runs')';
scraping_simpson_indices_runs = cell2mat(scraping_simpson_indices_runs')';

biopsy_shannon_indices_runs = cell2mat(biopsy_shannon_indices_runs')';
biopsy_simpson_indices_runs = cell2mat(biopsy_simpson_indices_runs')';

for i=1:N_runs % removes extra points, that are sometimes present if the simulation goes on for longer than 600 time units.
    if (length(biopsy_moranI_runs{i}) > 600)
        temp = biopsy_moranI_runs{i} ;
        temp = temp(1:600);
        biopsy_moranI_runs{i} = temp;
    end
end

biopsy_moranI_runs = cell2mat(biopsy_moranI_runs')';

for i=1:N_runs
    if (length(biopsy_gearyC_runs{i}) > 600)
        temp = biopsy_gearyC_runs{i} ;
        temp = temp(1:600);
        biopsy_gearyC_runs{i} = temp;
    end
end
                              
biopsy_gearyC_runs = cell2mat(biopsy_gearyC_runs')';

fpi0_runs = cell2mat(fpi0_vals')';
fpi_runs = cell2mat(fpi_vals')';
fpi2_runs = cell2mat(fpi2_vals')';
fpi3_runs = cell2mat(fpi3_vals')';

mpi_pos_runs = cell2mat(mpi_pos_vals')';
mpi_tot_runs = cell2mat(mpi_tot_vals')';

fpi0_bx_runs = cell2mat(fpi0_bx_vals')';
fpi_bx_runs = cell2mat(fpi_bx_vals')';
fpi2_bx_runs = cell2mat(fpi2_bx_vals')';
fpi3_bx_runs = cell2mat(fpi3_bx_vals')';
mpi_pos_bx_runs = cell2mat(mpi_pos_bx_vals')';
mpi_tot_bx_runs = cell2mat(mpi_tot_bx_vals')';
times_vals=cell2mat(time_of_cancer_occurrence_measure2_runs);

lattice_erpos_vals = cell2mat(lattice_erpos_runs')';
lattice_ki67_vals = cell2mat(lattice_ki67_runs')';

biopsy_erpos_vals = cell2mat(biopsy_erpos_runs')';
biopsy_ki67_vals = cell2mat(biopsy_ki67_runs')';

scraping_erpos_vals = cell2mat(scraping_erpos_runs')';
scraping_ki67_vals = cell2mat(scraping_ki67_runs')';

minTime = int64(fix(min(cell2mat(time_of_cancer_occurrence_measure2_runs)))); % This is the minimum time at which the first simulation reaches cancer


%% The following is the correction of the data for when we consider multiple bx at same time point
%% Specificially, it is to remove the 9 extra biopsies at the time point 50 --> only use this if using the code for multiple biopsies at a given timepoint
%% This code MUST be uncommented when using data for the case of multiple biopsies at a single timepoint and NOT plotting the multiple biopsy vs time data.
%                               
% biopsy_shannon_indices_runs = [biopsy_shannon_indices_runs(:,1:50) biopsy_shannon_indices_runs(:,60:end)];
% biopsy_simpson_indices_runs = [biopsy_simpson_indices_runs(:,1:50) biopsy_simpson_indices_runs(:,60:end)];
% biopsy_moranI_runs = [biopsy_moranI_runs(:,1:50) biopsy_moranI_runs(:,60:end)];
% biopsy_gearyC_runs = [biopsy_gearyC_runs(:,1:50) biopsy_gearyC_runs(:,60:end)];
% fpi0_bx_runs = [fpi0_bx_runs(:,1:50) fpi0_bx_runs(:,60:end)];
% fpi_bx_runs = [fpi_bx_runs(:,1:50) fpi_bx_runs(:,60:end)];
% fpi2_bx_runs = [fpi2_bx_runs(:,1:50) fpi2_bx_runs(:,60:end)];
% fpi3_bx_runs = [fpi3_bx_runs(:,1:50) fpi3_bx_runs(:,60:end)];
% biopsy_erpos_vals = [biopsy_erpos_vals(:,1:50) biopsy_erpos_vals(:,60:end)];
% biopsy_ki67_vals = [biopsy_ki67_vals(:,1:50) biopsy_ki67_vals(:,60:end)];
% 
% for i=1:N_runs
% temp  = biopsy_types_of_mutants_runs{i};
% temp = [temp(1:50);temp(60:end)];
%     biopsy_types_of_mutants_runs{i} = temp;
% end

                              
% The following is the code to make the LaTeX table of Hazards Ratios
varName = {'shannon_indices_runs','simpson_indices_runs','moranI_runs','gearyC_runs','fpi0_runs','fpi_runs','fpi2_runs','fpi3_runs','mpi_pos_runs','mpi_tot_runs','lattice_erpos_vals','lattice_ki67_vals',...
    'biopsy_shannon_indices_runs','biopsy_simpson_indices_runs','biopsy_moranI_runs','biopsy_gearyC_runs',...
    'fpi0_bx_runs','fpi_bx_runs','fpi2_bx_runs','fpi3_bx_runs','mpi_pos_bx_runs','mpi_tot_bx_runs','biopsy_erpos_vals','biopsy_ki67_vals'...
    'scraping_shannon_indices_runs','scraping_simpson_indices_runs','scraping_erpos_vals','scraping_ki67_vals' };
indexName ={ 'Shannon','Simpson','Moran''s I','Geary''s C','FPI 1','FPI 2',...
    'FPI 3','FPI 4','MPI Pos','MPI Tot','ER+','Ki-67','Shannon','Simpson','Moran''s I','Geary''s C','FPI 1','FPI 2',...
    'FPI 3','FPI 4','MPI Pos','MPI Tot','ER+','Ki-67','Shannon','Simpson','ER+','Ki-67'};
%for i=50:50:100
out = {};
for k=1:length(varName)
    [b1(k), ~,~,temp1] = eval(strcat('coxphfit(',varName{k},'(:,225),times_vals);')); % TIME SHOULD BE SET HERE FOR THE FIRST HR
    [b2(k), ~,~,temp2] = eval(strcat('coxphfit(',varName{k},'(:,225),times_vals);')); % TIME SHOULD BE SET HERE FOR THE SECOND HR
    
    r1(k) = temp1.p;
    r2(k) = temp2.p;
    z_score =1.96;
    unit_change = 0.1;
    conf_lower1 = exp(b1(k)  - (z_score*temp1.se));
    conf_lower2 = exp(b2(k)  - (z_score*temp2.se));
    conf_upper1 = exp(b1(k)  + (z_score*temp1.se));
    conf_upper2 = exp(b2(k)  + (z_score*temp2.se));
    conf_lower1 = conf_lower1 ^ unit_change;
    conf_lower2 = conf_lower2 ^ unit_change
    conf_upper1 = conf_upper1 ^ unit_change;
    conf_upper2 = conf_upper2 ^ unit_change
    
    prec_p =  2;
    prec_np = 2;
    %    p(k) = temp.p;
    b1(k) = exp(unit_change .* b1(k));
    b2(k) = exp(unit_change .* b2(k));
    if (r1(k) < 0.0001)
        p_val_1 = ' $< 10^{-4}$';
    else
        p_val_1 = strcat(' $ ',num2str(r1(k),prec_p),'$');
    end
    if (r2(k) < 0.0001)
        p_val_2 = ' $< 10^{-4}$';
        
    else
        p_val_2 = strcat(' $ ',num2str(r2(k),prec_p),'$');
        
    end
        
    p_val_1 = strcat(' \t &  $(' , num2str(conf_lower1,prec_p), ',' , num2str(conf_upper1,prec_p), ')$ \t & ',p_val_1);
    
    p_val_2 = strcat(' \t  &  $(' , num2str(conf_lower2,prec_p), ',' , num2str(conf_upper2,prec_p), ')$ \t & ',p_val_2);
    if (r1(k) < 0.05)
        first_str = strcat(' { \\bf \t ', num2str(b1(k),prec_np),' } ' );  

    else
        first_str =   num2str(b1(k),prec_np);
    end
    
    if (r2(k) < 0.05)
        second_str = strcat(' { \\bf \t ', num2str(b2(k),prec_np),' } ' );  

    else
        second_str =   num2str(b2(k),prec_np);
    end
    out{k} = strcat(' &  ',indexName{k},'  \t &   ',first_str,p_val_1, ' \t &  ' , second_str, p_val_2,'  \\\\ \n ');
    
end
outstr = '';
for i=1:length(out)
    outstr = strcat(outstr,' ',out{i});
end
sprintf(outstr)

%%--------- the following is the code to make the KM curves --------------------

% 
% r = figure();
% makeKMgraph(times_vals, shannon_indices_runs,80,'Shannon index','KM_lattice_shannon', r);
% 
% r = figure();
% makeKMgraph(times_vals, simpson_indices_runs,80,'Simpson index','KM_lattice_simpson', r);
% %set(gca,'title','Lattice Simpson Index')
% 
% r=  figure();
% makeKMgraph(times_vals, moranI_runs,80,'Moran''s \it I','KM_lattice_moranI', r);
% 
% r=  figure();
% makeKMgraph(times_vals, gearyC_runs,80,'Geary''s \it C','KM_lattice_gearyC', r);
% 
% r=  figure();
% makeKMgraph(times_vals, fpi0_runs,80,'FPI 1','KM_lattice_fpi0', r);
% 
% r=  figure();
% makeKMgraph(times_vals, fpi_runs,80,'PPI','KM_lattice_fpi', r);
% 
% r=  figure();
% makeKMgraph(times_vals, fpi2_runs,80,'FPI 3','KM_lattice_fpi2', r);
% 
% r=  figure();
% makeKMgraph(times_vals, fpi3_runs,80,'NPI','KM_lattice_fpi3', r);
% 
% %r=  figure();
% %makeKMgraph(times_vals, lattice_erpos_vals,80,'Proportion with at least 2 positive mutations','KM_lattice_erpos', r);
% 
% r=  figure();
% makeKMgraph(times_vals, lattice_ki67_vals,80,'Ki-67+ proportion','KM_lattice_ki67', r);
% 
% r = figure();
% makeKMgraph(times_vals, biopsy_shannon_indices_runs,80,'Shannon index','KM_bx_shannon', r);
% %set(h,'title','Lattice Shannon Index')
% 
% r = figure();
% makeKMgraph(times_vals, biopsy_simpson_indices_runs,80,'Simpson index','KM_bx_simpson', r);
% %set(gca,'title','Lattice Simpson Index')
% 
% r=  figure();
% makeKMgraph(times_vals, biopsy_moranI_runs,80,'Moran''s \it I','KM_bx_moranI', r);
% 
% r=  figure();
% makeKMgraph(times_vals, biopsy_gearyC_runs,80,'Geary''s \it C','KM_bx_gearyC', r);
% 
% r=  figure();
% makeKMgraph(times_vals, fpi0_bx_runs,80,'FPI 1','KM_bx_fpi0', r);
% 
% r=  figure();
% makeKMgraph(times_vals, fpi_bx_runs,80,'PPI','KM_bx_fpi', r);
% 
% r=  figure();
% makeKMgraph(times_vals, fpi2_bx_runs,80,'FPI 3','KM_bx_fpi2', r);
% 
% r=  figure();
% makeKMgraph(times_vals, fpi3_bx_runs,80,'NPI','KM_bx_fpi3', r);
% 
% %r=  figure();
% %makeKMgraph(times_vals, biopsy_erpos_vals,80,'Proportion with at least 2 positive mutations','KM_bx_erpos', r);
% 
% r=  figure();
% makeKMgraph(times_vals, biopsy_ki67_vals,80,'Ki-67+ proportion','KM_bx_ki67', r);
% 
% r = figure();
% makeKMgraph(times_vals, scraping_shannon_indices_runs,80,'Shannon index','KM_scraping_shannon', r);
% %set(h,'title','Lattice Shannon Index')
% 
% r = figure();
% makeKMgraph(times_vals, scraping_simpson_indices_runs,80,'Simpson index','KM_scraping_simpson', r);
% 
% %r=  figure();
% %makeKMgraph(times_vals, scraping_erpos_vals,80,'Proportion with at least 2 positive mutations','KM_scraping_erpos', r);
% 
% r=  figure();
% makeKMgraph(times_vals, scraping_ki67_vals,80,'Ki-67+ proportion','KM_scraping_ki67', r);

%-------------------------END of KM curve plotting section-----------------------------
                              
% -------- Correlation coefficient vs time graph plotting -----------------------------
% r_vals = zeros((minTime-1),10) ;
% r_vals_suppl_plot = zeros((minTime-1), 5); % just for ppi1, ppi , ppi3, mpipos, mpitot
%  for i=1:(minTime-1)
%     temp = corrcoef(shannon_indices_runs(:,i), times_vals);
%     r_vals(i,1) =temp(2,1);
% 
%     temp = corrcoef(simpson_indices_runs(:,i), times_vals);
%     r_vals(i,2) =temp(2,1);
% 
%     temp = corrcoef(moranI_runs(:,i), times_vals);
%     r_vals(i,3) =temp(2,1);
% 
%     temp = corrcoef(gearyC_runs(:,i), times_vals);
%     r_vals(i,4) =temp(2,1);
% 
%     temp = corrcoef(mpi_pos_runs(:,i), times_vals);
%     r_vals(i,5) =temp(2,1);
% 
%     temp = corrcoef(fpi_runs(:,i), times_vals);
%     r_vals(i,6) =temp(2,1);
% 
%     temp = corrcoef(mpi_tot_runs(:,i), times_vals);
%     r_vals(i,7) =temp(2,1);
% 
%     temp = corrcoef(fpi3_runs(:,i), times_vals);
%     r_vals(i,8) =temp(2,1);
% 
%     temp = corrcoef(lattice_erpos_vals(:,i), times_vals);
%     r_vals(i,9) =temp(2,1);
% 
%     temp = corrcoef(lattice_ki67_vals(:,i), times_vals);
%     r_vals(i,10) =temp(2,1);

% temp = corrcoef(fpi0_runs(:,i), times_vals);
% r_vals_suppl_plot (i,1) = temp(2,1);
% temp = corrcoef(fpi_runs(:,i), times_vals);
% r_vals_suppl_plot (i,2) = temp(2,1);
% temp = corrcoef(fpi2_runs(:,i), times_vals);
% r_vals_suppl_plot (i,3) = temp(2,1);
% temp = corrcoef(mpi_pos_runs(:,i), times_vals);
% r_vals_suppl_plot (i,4) = temp(2,1);
% temp = corrcoef(mpi_tot_runs(:,i), times_vals);
% r_vals_suppl_plot (i,5) = temp(2,1);
% end
% 
% 
% r = figure();
% r_vals_suppl_plot  = abs(r_vals_suppl_plot);
% plot(1:(minTime-1), r_vals_suppl_plot(:,1),'r','linewidth',2);
%  hold on
% plot(1:(minTime-1), r_vals_suppl_plot(:,2),'g','linewidth',2);
% plot(1:(minTime-1), r_vals_suppl_plot(:,3),'b','linewidth',2);
% plot(1:(minTime-1), r_vals_suppl_plot(:,4),'c','linewidth',2);
% plot(1:(minTime-1), r_vals_suppl_plot(:,5),'k','linewidth',2);
% legend('PPI 1', 'PPI', 'PPI 3', 'MPI 1' ,'MPI 2')
% r = figure();
%r_vals = abs(r_vals);
% plot(1:(minTime-1), r_vals(:,1),'r','linewidth',2);
% hold on
% plot(1:(minTime-1), r_vals(:,2),'g','linewidth',2);
% plot(1:(minTime-1), r_vals(:,3),'b','linewidth',2);
% plot(1:(minTime-1), r_vals(:,4),'k','linewidth',2);
% %plot(1:(minTime-1), r_vals(:,5),'r:','linewidth',2);
% plot(1:(minTime-1), r_vals(:,6),'g:','linewidth',2);
% %plot(1:(minTime-1), r_vals(:,7),'b:','linewidth',2);
% plot(1:(minTime-1), r_vals(:,8),'k:','linewidth',2);
% plot(1:(minTime-1), r_vals(:,9),'c','linewidth',2);
% plot(1:(minTime-1), r_vals(:,10),'m','linewidth',2);
% temp_val = legend('Shannon index','Simpson index','Moran''s \it I','Geary''s \it C', 'PPI','NPI','ER+ proportion','Ki-67+ proportion');%'FPI 1','FPI 2',' FPI 3', 'FPI 4','ER+ proportion' ,'Ki-67+ proportion');
% %xlabel ('Time','FontSize',24)
 %title ('Lattice','FontSize',24)
% 
 %xlim([1,minTime-1])
% %ylabel('Correlation coefficient','FontSize',24)
% set(gca,'FontSize',24)
%set(temp_val,'Location','SouthEastOutside')
% set(r,'Position',[100,100,1000,500])
% fName = 'corrCoeff_Lattice';
% print(r,'-depsc',[fName])
% savefig(r,fName)
% 
% -------- Code for Biopsy data graph----
                               
% r_vals = zeros((minTime-1),10) ;
% r_vals_suppl_plot = zeros((minTime-1), 5); % just for ppi1, ppi , ppi3, mpipos, mpitot

% for i=1:(minTime-1)
%     temp = corrcoef(biopsy_shannon_indices_runs(:,i), times_vals);
%     r_vals(i,1) =temp(2,1);
% 
%     temp = corrcoef(biopsy_simpson_indices_runs(:,i), times_vals);
%     r_vals(i,2) =temp(2,1);
% 
%     temp = corrcoef(biopsy_moranI_runs(:,i), times_vals);
%     r_vals(i,3) =temp(2,1);
% 
%     temp = corrcoef(biopsy_gearyC_runs(:,i), times_vals);
%     r_vals(i,4) =temp(2,1);
% 
%     temp = corrcoef(fpi0_bx_runs(:,i), times_vals);
%     r_vals(i,5) =temp(2,1);
% 
%     temp = corrcoef(fpi_bx_runs(:,i), times_vals);
%     r_vals(i,6) =temp(2,1);
% 
%     temp = corrcoef(fpi2_bx_runs(:,i), times_vals);
%     r_vals(i,7) =temp(2,1);
% 
%     temp = corrcoef(fpi3_bx_runs(:,i), times_vals);
%     r_vals(i,8) =temp(2,1);
% 
%     temp = corrcoef(biopsy_erpos_vals(:,i), times_vals);
%     r_vals(i,9) =temp(2,1);
% 
%     temp = corrcoef(biopsy_ki67_vals(:,i), times_vals);
%     r_vals(i,10) =temp(2,1);
% temp = corrcoef(fpi0_bx_runs(:,i), times_vals);
% r_vals_suppl_plot (i,1) = temp(2,1);
% temp = corrcoef(fpi_bx_runs(:,i), times_vals);
% r_vals_suppl_plot (i,2) = temp(2,1);
% temp = corrcoef(fpi2_bx_runs(:,i), times_vals);
% r_vals_suppl_plot (i,3) = temp(2,1);
% temp = corrcoef(mpi_pos_bx_runs(:,i), times_vals);
% r_vals_suppl_plot (i,4) = temp(2,1);
% temp = corrcoef(mpi_tot_bx_runs(:,i), times_vals);
% r_vals_suppl_plot (i,5) = temp(2,1);
%  end
% 
% 
%  r = figure();
% r_vals_suppl_plot  = abs(r_vals_suppl_plot);
% plot(1:(minTime-1), r_vals_suppl_plot(:,1),'r','linewidth',2);
%  hold on
% plot(1:(minTime-1), r_vals_suppl_plot(:,2),'g','linewidth',2);
% plot(1:(minTime-1), r_vals_suppl_plot(:,3),'b','linewidth',2);
% plot(1:(minTime-1), r_vals_suppl_plot(:,4),'c','linewidth',2);
% plot(1:(minTime-1), r_vals_suppl_plot(:,5),'k','linewidth',2);
% legend('PPI 1', 'PPI', 'PPI 3', 'MPI 1' ,'MPI 2')
%  r = figure();
% r_vals = abs(r_vals);
% plot(1:(minTime-1), r_vals(:,1),'r','linewidth',2);
% hold on
% plot(1:(minTime-1), r_vals(:,2),'g','linewidth',2);
% plot(1:(minTime-1), r_vals(:,3),'b','linewidth',2);
% plot(1:(minTime-1), r_vals(:,4),'k','linewidth',2);
% %plot(1:(minTime-1), r_vals(:,5),'r:','linewidth',2);
% plot(1:(minTime-1), r_vals(:,6),'g:','linewidth',2);
% %plot(1:(minTime-1), r_vals(:,7),'b:','linewidth',2);
% plot(1:(minTime-1), r_vals(:,8),'k:','linewidth',2);
% plot(1:(minTime-1), r_vals(:,9),'c','linewidth',2);
% plot(1:(minTime-1), r_vals(:,10),'m','linewidth',2);
% 
% temp_val = legend('Shannon index','Simpson index','Moran''s \it I','Geary''s \it C', 'PPI','NPI','ER+ proportion','Ki-67+ proportion');%'FPI 1','FPI 2',' FPI 3', 'FPI 4','ER+ proportion' ,'Ki-67+ proportion');
% %xlabel ('Time','FontSize',24)
%  title ('Biopsy','FontSize',24)
%  xlim([1,minTime-1])
% ylabel('Absolute correlation coefficient','FontSize',24)
% xlabel('Time','FontSize',24)
%  set(gca,'FontSize',24)
% set(temp_val,'Location','SouthEastOutside')
% set(r,'Position',[100,100,1000,500])
% fName = 'corrCoeff_Biopsy';
% print(r,'-depsc',[fName])
% savefig(r,fName)
% %%-------- Code for scraping based correlation coefficient graph ---------
% r_vals = zeros((minTime-1),4) ;
% for i=1:(minTime-1)
%     temp = corrcoef(scraping_shannon_indices_runs(:,i), times_vals);
%     r_vals(i,1) =temp(2,1);
% 
%     temp = corrcoef(scraping_simpson_indices_runs(:,i), times_vals);
%     r_vals(i,2) =temp(2,1);
% 
% 
%     temp = corrcoef(scraping_erpos_vals(:,i), times_vals);
%     r_vals(i,3) =temp(2,1);
% 
%     temp = corrcoef(scraping_ki67_vals(:,i), times_vals);
%     r_vals(i,4) =temp(2,1);
% end

% r = figure();
% r_vals = abs(r_vals);
% plot(1:(minTime-1), r_vals(:,1),'r','linewidth',2);
% hold on
% plot(1:(minTime-1), r_vals(:,2),'g','linewidth',2);
% plot(1:(minTime-1), r_vals(:,3),'b','linewidth',2);
% plot(1:(minTime-1), r_vals(:,4),'k','linewidth',2);
% 
% temp_val = legend('Shannon index','Simpson index','ER+ proportion' ,'Ki-67+ proportion');
% xlabel ('Time','FontSize',24)
% title ('Scraping','FontSize',24)
% 
% xlim([1,minTime-1])
% %ylabel('Correlation coefficient','FontSize',24)
% set(gca,'FontSize',24)
% set(temp_val,'Location','SouthEastOutside')
% set(r,'Position',[100,100,1000,500])
% fName = 'corrCoeff_Scraping';
% print(r,'-depsc',[fName])
% savefig(r,fName)

% -------- END of correlation coefficient vs time graph plotting -----------------------

% -------- Creating 2D plots of difference index vs later index ------------------------
plot_times_vec = 1:(minTime-1);
diff_corr_coeff = zeros(minTime-1,minTime-1,10);
diff_corr_coeff_2 = zeros(minTime-1,minTime-1,10);
diff_corr_coeff_3 = zeros(minTime-1,minTime-1,10);
varName = {'biopsy_shannon_indices_runs','biopsy_simpson_indices_runs','biopsy_moranI_runs','biopsy_gearyC_runs',...
    'fpi0_bx_runs','fpi_bx_runs','fpi2_bx_runs','fpi3_bx_runs','biopsy_erpos_vals','biopsy_ki67_vals'};
indexName  = {'Shannon index','Simpson index','Moran''s \it I','Geary''s \it C','FPI 1','IPP','FPI 3','INP','Proportion with at least 2 positive mutations','Mitotic proportion'};
for i=1:minTime-1
    for j=1:minTime-1
        if (i~=j)


            timeDiff = double(abs(i-j));

            for k = 1:length(varName)

                temp = eval(strcat('corrcoef ((', varName{k} , '(:,i) - ' , varName{k} , '(:,j)),times_vals);'));
                temp2 =eval(strcat('corrcoef(',varName{k} , '(:,min(i,j)),times_vals);')); % eval(strcat('corrcoef(',varName{k} , '(:,i),times_vals);')); %eval(strcat('corrcoef(',varName{k} , '(:,max(i,j)),times_vals);'));
                %    temp3 = eval(strcat('corrcoef ((',varName{k},'(:,i) - ',varName{k},'(:,j))./timeDiff,times_vals);'));
                diff_corr_coeff(i,j,k) = abs(temp(2,1));
                if (i<=j)
                diff_corr_coeff_2(i,j,k) = (abs(temp(2,1)) - abs(temp2(2,1))); %(abs(temp(2,1)) > abs(temp2(2,1))).* (abs(temp(2,1)) - abs(temp2(2,1)));
                else
                     diff_corr_coeff_2(i,j,k)=NaN;
                end
                %      diff_corr_coeff_3(i,j,k) = (abs(temp3(2,1)) > abs(temp2(2,1))).* (abs(temp3(2,1)) - abs(temp2(2,1)));
            end
            %

            %diff_corr_coeff_3(i,j) = (abss(temp3(2,1)) > abs(temp2(2,1))).* (abs(temp3(2,1)) - abs(temp2(2,1)));
        end

    end
end
fName_addon = {'shannon','simpson','moranI','gearyC','fpi0','fpi','fpi2','fpi3','erpos','ki67pos'};
for i=1:length(varName)
%     r= figure();
%     %subplot(10,2,i)
%     imagesc(diff_corr_coeff(:,:,i))
%     xlabel('Time 1','FontSize',24)
%     ylabel('Time 2','FontSize',24)
%     title(strcat('\Delta ',indexName{i}),'FontSize',24)
%     colorbar
%     set(gca,'FontSize',24)
%     colormap(jet)
%     fName = strcat('corr_coeff_bx_diff_',fName_addon{i});
%     print(r,'-depsc',[fName])
%     savefig(r,fName)


    r=figure();
    %subplot(10,2,i+1)
    imagesc(diff_corr_coeff_2(:,:,i))
    if (i==6)
    xlabel('Time 1','FontSize',24)
    ylabel('Time 2','FontSize',24)
        colorbar

    end
    title(strcat('\Delta ',indexName{i}),'FontSize',24)
    caxis([-0.4,0.15]);
    colormap(jet)
    set(gca,'FontSize',20)
set(gca,'YDir','normal');

    fName = strcat('diff_in_corr_coeff_bx_diff_',fName_addon{i});
    print(r,'-depsc',[fName])
    savefig(r,fName)

    %
    % r=figure();
    % imagesc(diff_corr_coeff_3(:,:,3))
    % xlabel('t_1','FontSize',20)
    % ylabel('t_2','FontSize',20)
    % title(strcat('Difference in correlation coefficient, \Delta ',indexName{i},' / \Delta t'),'FontSize',20)
    % colorbar
    % colormap(jet)
    %
    % print(figHandle,'-depsc',[fName])
    % savefig(figHandle,fName)


end
% -------- End of 2D plots of difference index vs later index ------------------------
                            
% ----------------Start of multiple biopsy at same timepoint graphs--------------------------------
%% Here, we take the 10 biopsy points taken at the time 50 and analyse them: we take out those biopsy indices and come up
% with the combi_avg, comb_max, etc over the 200 runs, and then see how that correlates with the end times

% varNames = {'biopsy_shannon_indices_runs','biopsy_simpson_indices_runs','biopsy_moranI_runs','biopsy_gearyC_runs',...
%     'fpi0_bx_runs','fpi_bx_runs','fpi2_bx_runs','fpi3_bx_runs','biopsy_erpos_vals','biopsy_ki67_vals'};

% indexName ={ 'Shannon index','Simpson index','Moran''s \it I','Geary''s \it C','FPI 1','PPI',...
%     'FPI 3','NPI','Proportion with at least 2 positive mutations','Mitotic proportion'};
% fName_addon = {'shannon','simpson','moranI','gearyC','fpi0','fpi','fpi2','fpi3','erpos','ki67pos'};

% numIndices=  length(varNames);
% num_rep_bx =10;
% ind_1 = 50;
% ind_2 = ind_1+ num_rep_bx-1;
% for j=1:numIndices
%     eval(strcat('index_j = ',varNames{j},' ;'));
%     index_j = index_j(:,ind_1:ind_2);
%     corrCoef_avg = [];
%     corrCoef_max = [];
%     corrCoef_min = [];
%     corrCoef_max_minus_min = [];
%     corrCoef_var = [];
% for i=1:num_rep_bx
%    index_j_avg = mean(index_j(:,1:i),2)
%    index_j_max = max(index_j(:,1:i),[],2);
%    index_j_min = min(index_j(:,1:i),[],2);
%    index_j_max_minus_min = index_j_max - index_j_min;
%   index_j_var  =  var(index_j(:,1:i),0,2);
%   temp = corrcoef(index_j_avg, times_vals);
%       corrCoef_avg (i) = temp(2,1);
%       temp = corrcoef(index_j_max, times_vals);
%       corrCoef_max (i) = temp(2,1);
%       temp = corrcoef(index_j_min, times_vals);
%       corrCoef_min (i) = temp(2,1);
%       temp =  corrcoef(index_j_max_minus_min, times_vals);
%       corrCoef_max_minus_min (i) = temp(2,1);
%       temp = corrcoef(index_j_var, times_vals);
%       corrCoef_var (i) = temp(2,1);
% end
% figHandle = figure();
% plot(1:num_rep_bx,abs(corrCoef_avg),'r','LineWidth',2)
% hold on
% plot(1:num_rep_bx,abs(corrCoef_max),'g','LineWidth',2)
% plot(1:num_rep_bx,abs(corrCoef_min),'b','LineWidth',2)
% plot(1:num_rep_bx,abs(corrCoef_max_minus_min),'k','LineWidth',2)
% plot(1:num_rep_bx,abs(corrCoef_var),'c','LineWidth',2)
%
% ylabel('Absolute Correlation coefficient','FontSize',24)
% xlabel('Number of biopsies','FontSize',24)
% xlim([1,10])
%
%  title(indexName{j}, 'FontSize',24)
%  fName =strcat('diffNumBx_ ',fName_addon{j});
%  set(gca, 'FontSize', 24)
%   print(figHandle,'-depsc',[fName])
%   savefig(figHandle,fName)
% end
% ----------------End of multiple biopsy at same timepoint graphs-------------------------------