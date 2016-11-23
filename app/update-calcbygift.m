%%% using gift.
%%% szu
%%% 2016-11-06

%%% set parameters
c2s_filenm = 'com2sub.txt';
%% set is881 as zero if 882 columns, first is compound index.
is881 = 1;
s2d_filenm = 'sub2domain.txt';
%% is2d_txt as zero if it's mat format, and save it as lamda_SubOnly
is2d_txt = 1;
%% file name for output, use .mat as the last chars for matlab format.
d2p_savefn = 'drug2protein.mat';

%%--------------------------------%%
%%% run, keep the following unchanged.
load('SUDO_Whole_Bio2012.mat');
load('Unique_Sub_Domain_140215.mat');
com2subs = importdata(c2s_filenm);
if is881 < 1
  com2subs = com2subs(:,2:end);
end
Drug2Sub = com2subs(:,SubIndex);
Drug2Sub = com2subs(:,UniqueColumn_sub);

if is2d_txt > 0
  Sub_Domain_Result = importdata(s2d_filenm);
  lamda_SubOnly = vec2mat(Sub_Domain_Result, length(UniqueColumn_dom));
  clear Sub_Domain_Result;
else
  load('s2d_filenm');
end
%% recover domain from unique to all.
Sub2Domain_Recover = zeros(size(lamda_SubOnly,1),length(size(Protein_Domain,2)));
for i = 1:length(ColumnContent_dom)
  for j = 1:length(ColumnContent_dom{i})
    Sub2Domain_Recover(:, ColumnContent_dom{i}(j)) = lamda_SubOnly(:,i);
  end
end
Sub2Domain_Recover(Sub2Domain_Recover>0.99) = 0.99;
%%-- Predict the Compound-Protein Interactions.
nrow = size(Drug2Sub,1);
ncol = size(Protein_Domain,1);
com2pro = zeros(nrow,ncol);
for i = 1:nrow
  subs = Drug2Sub(i,:);
  for j = 1:ncol
    doms = Protein_Domain(j,:);
    com2pro(i,j) = 1 - exp(sum(sum(log(1 - Sub2Domain_Recover(subs==1,doms==1)))));
  end
end
save(d2p_savefn,'com2pro');
