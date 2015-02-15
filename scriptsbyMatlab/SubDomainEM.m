function  [deltalamda,ZD] = SubDomainEM(maxstep,fn,fp,Drug_Sub,Drug_Protein,Protein_Domain,ZD) 

% GIFT by EM algorithm.
% Songpeng Zu / zusongpeng@gmail.com /
% 140226.
%% Load Data and variable declarations
[DrugNum,SubNum] = size(Drug_Sub);
[ProNum,DomNum] = size(Protein_Domain);
if nargin <7
    ZD = rand(SubNum,DomNum); % the probabilities of Sub and Domain Interact;
end
%% EM Algorithm

step = 1;
deltalamda = zeros(maxstep,1);
fprintf('%f\n',deltalamda(step));
while step <=maxstep % Using Delta ZD Matrix to be the standard
    tic
    step = step + 1
    lamda = ZD;
    % ---------------------------Calculate Expectation--------------------;
    % Calculate O;
    O = EStep(DrugNum,ProNum,fn,fp,Drug_Sub,Protein_Domain,lamda);
    %O = EStep_mex(DrugNum,ProNum,fn,fp,Drug_Sub,Protein_Domain,lamda);
    % -------------------------Calculate Maximum--------------------------;
    ZD = Mstep(SubNum,DomNum,Drug_Sub,Protein_Domain,Drug_Protein,lamda,fn,O);
    % ----------------------Judgement whether to be out of loop------------;
    deltalamda(step) = max(max(abs(ZD-lamda)))
    fprintf('%f\n',deltalamda(step));
    toc
end

end


