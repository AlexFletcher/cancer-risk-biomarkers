function makeROC_curve (index_vals, cancer_times)
% index_vals is a vector containing index values for the biomarker we are
% using
% cancer_times is a vector containing the time at which the simulation
% reaches cancer

N_sample_points = 250;
N_steps = 100;

%[b, ~,~,~] = coxphfit(index_vals,cancer_times);
%index_vals=exp(b.*index_vals);
max_index_val = max(index_vals');

min_time = min(cancer_times);
max_time = max(cancer_times);
time_cutoffs = min_time:(max_time-min_time)/N_steps: max_time;

cutoffs = 0:(max_index_val/N_sample_points):max_index_val;
ROC = zeros((length(cutoffs)),2,length(time_cutoffs));

%figure();

AUC = zeros(length(time_cutoffs),1);
ppv = zeros(length(cutoffs), length(time_cutoffs));
npv = zeros(length(cutoffs), length(time_cutoffs));

for j=1:length(time_cutoffs)
     cancer_times_cut = cancer_times <= time_cutoffs(j);


    for i=1:length(cutoffs)
        
        index_vals_cut = index_vals >= cutoffs(i);

        tp = sum(cancer_times_cut .* index_vals_cut);
        fn = sum(cancer_times_cut .* (1-index_vals_cut));
        tn = sum((1-cancer_times_cut) .* (1-index_vals_cut));
        fp = sum((1-cancer_times_cut) .* (index_vals_cut));
        sens = tp/(tp+fn);
        spec = tn/(tn+fp);
        ppv(i,j) = tp/(tp+fp);
        npv(i,j) = tn/(tn+fn);
        
        if (~isnan(sens) && (~isnan(spec)))
            ROC(i,:,j)  = [(1-spec) sens];
        else
            ROC(i,:,j)=[0 0];
        end
    end
    [x, y] = sort(ROC(:,1,j),1);
    ROC(:,:,j)= ROC(y,:,j);
   % plot(ROC(:,1,j), ROC(:,2,j),'-')
    AUC(j) = trapz(ROC(:,1,j), ROC(:,2,j));
    %hold on
end

[maxAUC, max_ind]= max(AUC)

x_vals=ROC(:,1,max_ind);
y_vals=ROC(:,2,max_ind);

max(ppv)
max(npv)

%plot(x_vals,y_vals,'r.');
plot(time_cutoffs, AUC)
hold on
end