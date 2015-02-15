function varSub2Domain = devEM(fn,fp,Sub2Domain,likehood,quesidev,Drug2Protein,Drug2Sub,Protein2Domain)

% Get the variance of Sub2Domain Matrix.
% Songpeng Zu
% 2014-07-13

%% Get the result.
likehood2 = likehood.^2;
varSub2Domain = zeros(size(Sub2Domain,1),size(Sub2Domain,2));
for i = 1:size(Sub2Domain,1)
    DrugIndex = find(Drug2Sub(:,i) == 1);
    for j = 1:size(Sub2Domain,2)
        ProIndex = find(Protein2Domain(:,j) == 1);
        t1 = Drug2Protein(DrugIndex,ProIndex)./likehood2(DrugIndex,ProIndex);
        t2 = (1 - Drug2Protein(DrugIndex,ProIndex))./...
            (1-likehood(DrugIndex,ProIndex)).^2;
        t = t1(:)+t2(:);
        %s = (1-fn-fp)*exp(sum(sum(quesidev(DrugIndex,ProIndex)-log(1-Sub2Domian(i,j))))); % Note: How to avoid Sub2Domain(i,j)==1
        s = (1-fn-fp)*quesidev(DrugIndex,ProIndex)/(1-Sub2Domain(i,j));
        varSub2Domain(i,j) = 1/dot(s(:).^2,t);
    end
end


end

