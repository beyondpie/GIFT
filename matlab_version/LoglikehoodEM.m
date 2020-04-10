function loglikehood = LoglikehoodEM(DrugNum,ProNum,Drug_Sub,Protein_Domain,ZD,Drug_Protein,fn,fp)

% This function is to calculate Loglikehood in EM.
% Songpeng Zu
% 140226

%% Main Function
loglikehood = 0;
for j=1:ProNum 
        ProIndex  = Protein_Domain(j,:)==1;
        for  i=1:DrugNum                                 
            YP = 1 - exp(sum(sum(log(1-ZD(Drug_Sub(i,:)==1,ProIndex)))));
            O_ij = (1-fn)*YP + fp*(1-YP);
             if Drug_Protein(i,j)==1
                 loglikehood = loglikehood + log10(O_ij);
             else
                 loglikehood = loglikehood + log10(1-O_ij);
             end
        end
end


end

