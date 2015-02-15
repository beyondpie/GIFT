function Sub2Domain = AssociationMethod(Drug2Protein,Protein2Domain,Drug2Sub)
%This function is to get the Sub2Domain Matrix by Association Method.
%   INPUT:Drug2Protein, Protein2Domain, Drug2Sub
%   OUTPUT:Sub2Domain
%   Songpeng Zu /*zusongpeng@gmail.com*/
%   2014-02-22

Sub2Domain = zeros(size(Drug2Sub,2),size(Protein2Domain,2));
for s = 1:size(Drug2Sub,2)
    DrugIndex = find(Drug2Sub(:,s)==1);
    for d = 1:size(Protein2Domain,2)    
        ProteinIndex = find(Protein2Domain(:,d)==1);
        I = length(find(Drug2Protein(DrugIndex,ProteinIndex)==1));
        Sub2Domain(s,d) = I/(length(DrugIndex)*length(ProteinIndex));
    end
end

end

