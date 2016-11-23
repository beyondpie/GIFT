function result = GIFTinloop(comi,proj)
    global Drug2Sub Protein_Domain Sub2Domain_Recover
    subs = Drug2Sub(comi,:);
    doms = Protein_Domain(proj,:);
    result = 1 - exp(sum(sum(log(1 - Sub2Domain_Recover(subs==1,doms==1)))));
end
