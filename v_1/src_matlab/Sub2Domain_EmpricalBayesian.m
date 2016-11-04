function Sub2Domain_EB_Init = Sub2Domain_EmpricalBayesian(Drug2Protein,Protein_Domain,Drug_Sub_Reduce)
%This function is to get the Sub2Domain_EM matrix based on
%Sub2Domain_Empirical Bayesian method.
%   INPUT:Drug2Protein,Protein2Domain,Drug2Sub
%   OUTPUT:Sub2Domain_EB
%   Songpeng Zu /*zusongpeng@gmail.com*/
%   2014-02-22

tlinenum = size(Drug_Sub_Reduce,2)*size(Protein_Domain,2);
Sub2Domain = zeros(tlinenum,5);
for s = 1:size(Drug_Sub_Reduce,2);
    for d = 1:size(Protein_Domain,2)
        currentlinenum = (s-1)*size(Protein_Domain,2) + d;
        Sub2Domain(currentlinenum,1) = s;
        Sub2Domain(currentlinenum,2) = d;
        sline = find(Drug_Sub_Reduce(:,s)==1);
        dline = find(Protein_Domain(:,d)==1);
        Sub2Domain(currentlinenum,3) = length(sline)*length(dline);
        I = length(find(Drug2Protein(sline,dline)==1));
        Sub2Domain(currentlinenum,4) = I;
        Sub2Domain(currentlinenum,5) = Sub2Domain(currentlinenum,4)/Sub2Domain(currentlinenum,3);
    end
end

% Theta <-Beta(alpha,beta), Estimate the alpah,beta with mean and variance
% of data.
cutoff = mean(sum(Drug_Sub_Reduce))*mean(sum(Protein_Domain));
data_tmp = Sub2Domain(Sub2Domain(:,3)>cutoff,:);
theta_mean = mean(data_tmp(:,5));
theta_var = var(data_tmp(:,5));

%theta_mean = mean(Sub2Domain(:,5));
%theta_var = var(Sub2Domain(:,5));
 plus_tmp = theta_mean*(1-theta_mean)/theta_var-1;
 
 alpha = theta_mean*plus_tmp;
 beta = plus_tmp - alpha;
 % Revise the AS Score by mean of posterior.
 Sub2Domain_EB = zeros(size(Sub2Domain,1),1);
 for i = 1:length(Sub2Domain_EB)
     Sub2Domain_EB(i) = (alpha+Sub2Domain(i,4))/...
         (alpha+beta+Sub2Domain(i,3));
 end
 Sub2Domain(:,6) = Sub2Domain_EB;
% Write to a New Lamda Inition Matrix;
Sub2Domain_EB_Init = zeros(size(Drug_Sub_Reduce,2),size(Protein_Domain,2));
for s = 1:size(Drug_Sub_Reduce,2)
    for d = 1:size(Protein_Domain,2)
        currentlinenum = (s-1)*size(Protein_Domain,2) + d;
        Sub2Domain_EB_Init(s,d) = Sub2Domain(currentlinenum,6);
    end
end


end

