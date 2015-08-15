function makeKMgraph(times, comp_vals, spec_time_index,title_var,fName,figHandle)

% This code creates the Kaplan-Meier survival curve, by separating the times variable by the quartiles of the comp_vals vector, when the index value is taken at spec_time_index
% Also, creates a csv with the file name specified by fName, for input into the R code to perform the logrank test.
% Note that in order for R to use this csv, the file must be manually modified to have the header "end_time" for column 1 and "quartile" for column 2

comp_vals = comp_vals(:,spec_time_index);

cutoffs = prctile(comp_vals,[25 50 75]);

nxt_uq = 1;
nxt_umq = 1;
nxt_lmq = 1;
nxt_lq = 1;

uq=[];
umq = [0];
lmq = [0];
lq = [0];

quartile = zeros(length(times) , 1);

for i=1:length(times)
   if (comp_vals(i) > cutoffs(3))
        uq(nxt_uq) = times(i);
        nxt_uq = nxt_uq+1;
        quartile(i) = 1;
   elseif(comp_vals(i) > cutoffs(2))
        umq(nxt_umq) = times(i);
        nxt_umq = nxt_umq+1;
        quartile(i) = 2;

   elseif(comp_vals(i) > cutoffs(1))
        lmq(nxt_lmq) = times(i);
        nxt_lmq = nxt_lmq+1;
        quartile(i) = 3;
   else
       lq(nxt_lq) = times(i);
       nxt_lq = nxt_lq+1;
quartile(i) = 4;

   end       
    
end



ecdf(lq,'function','survivor')

hold on
ecdf(lmq,'function','survivor')
ecdf(umq,'function','survivor')
ecdf(uq,'function','survivor')

title(title_var,'FontSize',24);
set(gca,'FontSize',24) %'LineWidth',2)
h  = get(gca,'Children')
h = findobj('type','line')
set(h(1),'Color','r','LineWidth',2);
set(h(2),'Color','g','LineWidth',2);
set(h(3),'Color','b','LineWidth',2);
set(h(4),'Color','k','LineWidth',2);
xlabel(' ') %Time
ylabel(' ') %Probability of no cancer
set(figHandle,'Position',[100,100,500,400])

print(figHandle,'-depsc',[fName])
savefig(figHandle,fName)

csvMat = [times quartile]
csvwrite(strcat(fName,'.csv'),csvMat)


%--- the following is the code that is for the case of wanting to do upper
%half vs lower half, and not quartiles -----

% cutoffs = prctile(comp_vals,[50]);
% 
% nxt_uq = 1;
% nxt_umq = 1;
% nxt_lmq = 1;
% nxt_lq = 1;
% 
% uq=[];
% umq = [0];
% lmq = [0];
% lq = [0];
% quartile = zeros(length(times) , 1);
% for i=1:length(times)
%    if (comp_vals(i) >= cutoffs(1))
%        uq(nxt_uq) = times(i);
%        nxt_uq = nxt_uq+1;
% quartile(i) = 1;
% %    elseif(comp_vals(i) > cutoffs(1))
% %       umq(nxt_umq) = times(i);
% %       nxt_umq = nxt_umq+1;
% % quartile(i) = 2;
% % 
% %   elseif(comp_vals(i) > cutoffs(1))
% %       lmq(nxt_lmq) = times(i);
% %       nxt_lmq = nxt_lmq+1;
% % quartile(i) = 3;
% 
%    else
%        lq(nxt_lq) = times(i);
%        nxt_lq = nxt_lq+1;
% quartile(i) = 2;
% 
%    end       
%     
% end
% 
% 
% 
% ecdf(lq,'function','survivor')
% hold on
% ecdf(lmq,'function','survivor')
% ecdf(umq,'function','survivor')
% ecdf(uq,'function','survivor')

% title(title_var,'FontSize',24);
% set(gca,'FontSize',24) %'LineWidth',2)
% h  = get(gca,'Children')
% h = findobj('type','line')
% set(h(1),'Color','r','LineWidth',2);
% set(h(2),'Color','g','LineWidth',2);
% set(h(3),'Color','b','LineWidth',2);
% set(h(4),'Color','k','LineWidth',2);
% xlabel('Time')
% ylabel('Probability of no cancer')
% set(figHandle,'Position',[100,100,500,400])
% 
% print(figHandle,'-depsc',[fName])
% savefig(figHandle,fName)
% 
% csvMat = [times quartile]
% csvwrite(strcat(fName,'.csv'),csvMat)

end


