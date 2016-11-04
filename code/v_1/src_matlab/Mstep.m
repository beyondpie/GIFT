function ZD = Mstep(SubNum,DomNum,Drug_Sub,Protein_Domain,Drug_Protein,lamda,fn,O)%#codegen

%   This function is to establish Mstep.
%   Songpeng Zu
%   140226

%% Main Function
ZD = zeros(SubNum,DomNum);
for m=1:SubNum
        Y = sum(Drug_Sub(:,m)==1);
        DrugIndex = find(Drug_Sub(:,m)==1);
        for n=1:DomNum
             ProteinIndex = find(Protein_Domain(:,n)==1);
             DI = repmat(DrugIndex,length(ProteinIndex),1);
             PI = (ProteinIndex*ones(1,length(DrugIndex)))';
             PI = PI(:);
             OneI = find(Drug_Protein(sub2ind(size(Drug_Protein),DI,PI))==1);
             ZeroI = find(Drug_Protein(sub2ind(size(Drug_Protein),DI,PI))==0);
             temp = sum((1-fn)./O(sub2ind(size(O),DI(OneI),PI(OneI)))) + sum(fn./(1-O(sub2ind(size(O),DI(ZeroI),PI(ZeroI)))));
             ZD(m,n) = lamda(m,n)*temp/(Y*sum(Protein_Domain(:,n)==1));
             %clear ProteinIndex DI PI OneI ZeroI
        end
        %clear DrugIndex
    
end

