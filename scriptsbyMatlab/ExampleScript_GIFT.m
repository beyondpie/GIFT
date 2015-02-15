% Example code for GIFT.
% Songpeng Zu /*zusongpeng@gmail.com*/
% 2014-11-28

%% Read Original Data.
% The original data are downloaded from http://cbio.ensmp.fr/~yyamanishi/l1binary/ 
% Suppose you:
% save the drug representation matrix with PubChem chemical substructures as "drug_repmat.txt"
% save the target protein representation matrix with PFAM domains as "target_repmat.txt"
% save the adjacency matrix for drug-target interactions as "inter_admat.txt"

% Drug Sub Data
Drug_Sub_Struct = importdata('drug_repmat.txt');
Drug_Sub = Drug_Sub_Struct.data(:,2:end); % This is the drug2sub matrix 
DrugCID = Drug_Sub_Struct.data(:,1); 

% Protein Domain Data
Protein_Domain_Struct = importdata('target_repmat.txt');
Protein_Domain = Protein_Domain_Struct.data; % This is the protein2domain matrix 

% Drug Protein Data
Drug_Protein_Struct = importdata('inter_admat.txt');
Drug_Protein = Drug_Protein_Struct.data(:,2:end); % This is the drug2protein matrix

% You can save the three matrix above, or you can follow the basic format to organize your own data set.
% These data can be saved for the C++ version of GIFT.
% You can use the matlab function, or you can run the command like:
write_matrix_int_txt('filename',Drug_Sub);

%% Initialize the substructure_domain score matrix for GIFT.
Sub2Domain_EmpBayes = Sub2Domain_EmpricalBayesian(Drug_Protein,Protein_Domain,Drug_Sub);
% This the initial sub2domain matrix for GIFT. You can save it as textfile for the C++ version of GIFT
write_matrix_double_txt('filename', Sub2Domain_EmpBayes);

%% The Association Method can be used to get a simple and local estimation of Substructure-Domain Scores.
Sub2Domain_AS = AssociationMethod(Drug_Protein,Protein_Domain,Drug_Sub);

%% Running GIFT by Matlab
% Note: We write the matlab version of GIFT. But since this algorithm is an
% iterative algorithm, its speed is relatively slow. We suggest you
% prepare the data following the steps above and save them. Then use our
% C++ version of GIFT, which is much faster.

fn = 0.85; % the false negative rate 
fp = 0.0001; % the false positive rate
% In GIFT, we use maxstep number to control the stop instead of calculating
% the delta of loglikehood due to the speed. You can modify it. Based on our practive, 
% It converges well in current data set if you set it as 600. 
maxstep = 600;
% Then you can use the SubDomainEM function to run GIFT.
[deltalamda,ZD] = SubDomainEM(maxstep,fn,fp,Drug_Sub,Drug_Protein,Protein_Domain,Sub2Domain_EmpBayes);
% deltalamda record the max absolute delta between the matrix in the previous step and the one now.
% deltalamda is used to help us to control the convergence.
% ZD is the Final result by GIFT.

%% Running GIFT by GIFT of C++ version
% You just need to save the data above.
% Then aftering install the packages (gsl and boost) needed for GIFT, you
% can compile it. Then use the command line to run it.
% The details can be found in the user manual. 

% You might need some time to wait for it. Then, the result will be saved as one txt file.
% You can analyse it by Matlab.

% Read the result.
tmpS2D = importdata(resultfilepath); 
% The matrix is saved by gsl function following the column order.
% You can convert it as the real matrix by the function below.
S2D = vec2mat(tmpS2D,size(Protein_Domain,2));  % size(Protein_Domain,2) is the number of domains.

%% Variance estimation.
% We use the observed Fisher information to evalute the variances of the results.

% Before calculating the variances, we need to get the likelihood matrix and the quesi-derivation of drug_protein Matrix. 

% The likelihood matrix
[DrugNum,ProNum] = size(Drug_Protein);
O = zeros(DrugNum,ProNum);
for j = 1:ProNum
    ProIndex = Protein_Domain(j,:) == 1;
    for i = 1:DrugNum
        YP = 1 - exp(sum(sum(log(1-ZD(Drug_Sub(i,:)==1,ProIndex)))));
        O(i,j) = (1-fn)*YP + fp*(1-YP);
    end
end

% The quesi-derivation of drug_protein Matrix.
[DrugNum,ProNum] = size(Drug_Protein);
quesidev = zeros(DrugNum,ProNum);                                   
Sub2Domain = ZD;
Sub2Domain(ZD>=0.95) = 0.95;
for j = 1:ProNum
    ProIndex = Protein_Domain(j,:) == 1;
    for i = 1:DrugNum
        quesidev(i,j) = exp(sum(sum(log(1-Sub2Domain(Drug_Sub(i,:)==1,ProIndex)))));
    end
end   

% Then we can get the variance estiamtion using the devEM function.
varSub2Domain = devEM(fn,fp,Sub2Domain,O,quesidev,Drug_Protein,Drug_Sub,Protein_Domain);

%% Get the significant substructure-domain interactions.
% If the predicted score is larger than its standrad variance, we treat it as a significant substructure-domain interaction. 
