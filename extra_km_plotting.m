
% % %% this is the code to get the proportion of each cell type # of pos mutations
num_times = length(types_of_mutants_time_series_runs{1});
N_runs = length(types_of_mutants_time_series_runs);
N_size =100;
max_mut_at_time = zeros((minTime-1),1);
    for j=1:(minTime-1)
            max_val = -9999;

        for i=1:N_runs

        temp_val = types_of_mutants_time_series_runs{i};
        temp_val = temp_val{j};
        % now TemP_val is the types of mutatns matrix at time i

        cur_val = max(temp_val(:,3));
        if (cur_val > max_val)
            max_val= cur_val;
        end

    end
    max_mut_at_time(j) = max_val;
end

highest_max = max(max_mut_at_time(:));

props = zeros(minTime-1,highest_max);
for k=1:(minTime-1)

for i=1:N_runs
    for j=1:highest_max
            temp_val = types_of_mutants_time_series_runs{i};
            types_of_mutants = temp_val{k};

            props(k,j)=props(k,j) + (sum(types_of_mutants(:,4).*(types_of_mutants(:,3)== j))/(N_size*N_size));
        end
end
%            props(k,j) = props(k,j)./N_runs;

end

props = props ./ N_runs;

plot(1:(minTime-1),props(:,1),'r','LineWidth',2)
hold on
plot(1:(minTime-1),props(:,2),'g','LineWidth',2)
plot(1:(minTime-1),props(:,3),'b','LineWidth',2)
plot(1:(minTime-1),props(:,4),'k','LineWidth',2)
plot(1:(minTime-1),props(:,5),'c','LineWidth',2)

% plot(1:(minTime-1),props(:,6),'r:','LineWidth',2)
% plot(1:(minTime-1),props(:,7),'g:','LineWidth',2)
% plot(1:(minTime-1),props(:,8),'b:','LineWidth',2)
% plot(1:(minTime-1),props(:,9),'k:','LineWidth',2)
% plot(1:(minTime-1),props(:,10),'c:','LineWidth',2)
% 
% plot(1:(minTime-1),props(:,11),'r-.','LineWidth',2)
% plot(1:(minTime-1),props(:,12),'g-.','LineWidth',2)
% plot(1:(minTime-1),props(:,13),'b-.','LineWidth',2)
% plot(1:(minTime-1),props(:,14),'k-.','LineWidth',2)
% plot(1:(minTime-1),props(:,15),'c-.','LineWidth',2)

legend('1','2','3','4','5')%,'6','7','8','9','10','11','12','13','14','15')
xlim([1,(minTime-1)])
xlabel('Time','FontSize',24)
ylabel('Proportion of cells','FontSize',24)
set(gca,'FontSize',24)
title('Proportion of cells with given number of positive mutations','FontSize',24)

% % now should be able to plot the props vector and show the different lines
% % over time...

                            % %
% lattice_props_of_cell_muts = zeros(N_runs,minTime-1,N_m);
% lattice_indicator_of_cell_muts = zeros(N_runs,minTime-1,N_m);
% biopsy_props_of_cell_muts = zeros(N_runs,minTime-1,N_m);
% biopsy_indicator_of_cell_muts = zeros(N_runs,minTime-1,N_m);
% scraping_props_of_cell_muts = zeros(N_runs,minTime-1,N_m);
% scraping_indicator_of_cell_muts = zeros(N_runs,minTime-1,N_m);
% 
% corrCoeff_lattice_prop = zeros(minTime-1,N_m);
% corrCoeff_lattice_ind = zeros(minTime-1,N_m);
% corrCoeff_biopsy_prop = zeros(minTime-1,N_m);
% corrCoeff_biopsy_ind = zeros(minTime-1,N_m);
% corrCoeff_scraping_prop = zeros(minTime-1,N_m);
% corrCoeff_scraping_ind = zeros(minTime-1,N_m);
% 
% cutoff_prop = 0.05;
% N_size = 100;%% these vars should be moved to the top!
% N_bx_size = 20;
% N_scraping = 1000;
% for i=1:N_runs
%     for j = 1:(minTime-1)
%         for k= 1:N_m
%             temp_mat = types_of_mutants_time_series_runs{i};
%             temp_mat = temp_mat{j};
%             lattice_props_of_cell_muts(i,j,k) = (sum((temp_mat(:,3) >= k).*(temp_mat(:,4)))) ./ N_size.*N_size;
%             lattice_indicator_of_cell_muts(i,j,k) = lattice_props_of_cell_muts(i,j,k) >= cutoff_prop ;
%             temp_mat = biopsy_types_of_mutants_runs{i};
%             temp_mat = temp_mat{j};
%             biopsy_props_of_cell_muts(i,j,k) = (sum((temp_mat(:,3) >= k).*(temp_mat(:,4)))) ./ N_bx_size.*N_bx_size;
%             biopsy_indicator_of_cell_muts(i,j,k) = biopsy_props_of_cell_muts(i,j,k) >= cutoff_prop ;
%             temp_mat = scraping_types_of_mutants_runs{i};
%             temp_mat = temp_mat{j};
%             scraping_props_of_cell_muts(i,j,k) = (sum((temp_mat(:,3) >= k).*(temp_mat(:,4)))) ./ N_scraping;
%             scraping_indicator_of_cell_muts(i,j,k) = scraping_props_of_cell_muts(i,j,k) >= cutoff_prop ;
%         end
%     end
% end
% 
% for i=1:(minTime-1)
%     for k = 1:N_m
%         temp = corrcoef(lattice_props_of_cell_muts(:,i,k),times_vals);
%         corrCoeff_lattice_prop(i,k) = temp(2,1);
%         temp = corrcoef(lattice_indicator_of_cell_muts(:,i,k),times_vals);
%         corrCoeff_lattice_ind(i,k) = temp(2,1);
%         temp =  corrcoef(biopsy_props_of_cell_muts(:,i,k),times_vals);
%         corrCoeff_biopsy_prop(i,k) =temp(2,1);
%         temp = corrcoef(biopsy_indicator_of_cell_muts(:,i,k),times_vals);
%         corrCoeff_biopsy_ind(i,k) = temp(2,1);
%         temp =  corrcoef(scraping_props_of_cell_muts(:,i,k),times_vals);
%         corrCoeff_scraping_prop(i,k) =temp(2,1);
%         temp = corrcoef(scraping_indicator_of_cell_muts(:,i,k),times_vals);
%         corrCoeff_scraping_ind(i,k) = temp(2,1);
%     end
% end
% for k=1:N_m
%     plot(1:(minTime-1),corrCoeff_lattice_ind(:,k),'color',rand(1,3));
%     hold on
% end
% legend('1','2','3','4','5','6','7','8','9','10')
% 
% 
% 
% %% the following is the code to make the latex table
% varName = {'lattice_props_of_cell_muts','biopsy_props_of_cell_muts','scraping_props_of_cell_muts','lattice_indicator_of_cell_muts','biopsy_indicator_of_cell_muts','scraping_indicator_of_cell_muts'};
% %indexName ={ 'Lattice, proportions','Biopsy, proportions','Scraping, proportions'};
% %for i=50:50:100
% out = {};
% for k=1:length(varName)
%     for j=1:N_m
%         
%         [b1(k), ~,~,temp1] = eval(strcat('coxphfit(',varName{k},'(:,50,', num2str(j),'),times_vals);'));
%         [b2(k), ~,~,temp2] = eval(strcat('coxphfit(',varName{k},'(:,100,', num2str(j),'),times_vals);'));
%         
%         r1(k) = temp1.p;
%         r2(k) = temp2.p;
%         z_score =1.96;
%         unit_change = 1;
%         conf_lower1 = exp(b1(k)  - (z_score*temp1.se));
%         conf_lower2 = exp(b2(k)  - (z_score*temp2.se));
%         conf_upper1 = exp(b1(k)  + (z_score*temp1.se));
%         conf_upper2 = exp(b2(k)  + (z_score*temp2.se));
%         conf_lower1 = conf_lower1 ^ unit_change;
%         conf_lower2 = conf_lower2 ^ unit_change;
%         conf_upper1 = conf_upper1 ^ unit_change;
%         conf_upper2 = conf_upper2 ^ unit_change;
%         
%         prec_p =  2;
%         prec_np = 2;
%         %    p(k) = temp.p;
%         b1(k) = exp(unit_change .* b1(k));
%         b2(k) = exp(unit_change .* b2(k));
%         if (r1(k) < 0.0001)
%             p_val_1 = ' $< 10^{-4}$';
%         else
%             p_val_1 = strcat(' $ ',num2str(r1(k),prec_p),'$');
%         end
%         if (r2(k) < 0.0001)
%             p_val_2 = ' $< 10^{-4}$';
%             
%         else
%             p_val_2 = strcat(' $ ',num2str(r2(k),prec_p),'$');
%             
%         end
%         p_val_1 = strcat(' &  $(' , num2str(conf_lower1,prec_p), ',' , num2str(conf_upper1,prec_p), ')$ & ',p_val_1);
%         
%         p_val_2 = strcat(' &  $(' , num2str(conf_lower2,prec_p), ',' , num2str(conf_upper2,prec_p), ')$ & ',p_val_2);
%         
%         out{j,k} = strcat(' &  ',num2str(j),' Mutations   &   ',num2str(b1(k),prec_np),p_val_1, '  &  ' , num2str(b2(k),prec_np), p_val_2,'  \\\\ \n ');
%         
%     end
% end
% 
%     outstr = '';
% 
%     for i=1:length(out(1,:))
%         for j=1:length(out(:,1))
%         outstr = strcat(outstr,' ',out{j,i});
%         end
%     end
%     sprintf(outstr)
