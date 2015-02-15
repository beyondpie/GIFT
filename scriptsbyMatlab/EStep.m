function O = EStep(DrugNum,ProNum,fn,fp,Drug_Sub,Protein_Domain,lamda) %#codegen

% This function is to establish Estep.
% Songpeng Zu
% 140226

%% Main Function

O = zeros(DrugNum,ProNum);
for j = 1:ProNum
    ProIndex = Protein_Domain(j,:) == 1;
    for i = 1:DrugNum
        YP = 1 - exp(sum(sum(log(1-lamda(Drug_Sub(i,:)==1,ProIndex)))));
        O(i,j) = (1-fn)*YP + fp*(1-YP);
    end
end

end

