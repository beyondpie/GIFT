% Calculate the targets of 800 compounds from TCM by GIFT.
% Songpeng Zu
% 2015-04-09

%% Load data.
% tcm2subs = textread('tcm2subs800_pubchem.txt'); % 882 columns
% drug2subs = textread('ANdrugs2sub.txt');
% com2subs = [tcm2subs(:,2:end);drug2subs(:,2:end)];

com2subs = importdata('SubStructure_111110_Final.txt'); % 881 columns


%coms = [tcm2subs(:,1);drug2subs(:,1)];
%write_matrix_int_txt('compoundslistcid.txt',coms);
%% Predicte the targets.
%-- Get the Sub2Domain Matrix first.
% Here we use SubOnly without the combinations of subs.
load('SUDO_Whole_Bio2012.mat');
%load('SubCom_660_140222.mat');
load('Unique_Sub_Domain_140215.mat');
proteinnames = textread('target_UniID','%s');

DomNum = length(UniqueColumn_dom);
Sub_Domain_Result = importdata('ReduceSubOnly\Sub_Domain_Result.txt');
lamda_SubOnly = vec2mat(Sub_Domain_Result,DomNum);
clear Sub_Domain_Result;

global Drug2Sub Protein2Domain Sub2Domain_Recover
Drug2Sub = com2subs(:,SubIndex);
Drug2Sub = com2subs(:,UniqueColumn_sub);

Protein2Domain = Protein_Domain(:,UniqueColumn_dom);
domPF = textread('domPF.txt','%s');

Sub2Domain_Recover = zeros(size(lamda_SubOnly,1),length(domPF));
## for i = 1:length(ColumnContent_dom)
##     Sub2Domain_Recover(:,ColumnContent_dom{i}(1)) = lamda_SubOnly(:,i);
##     if length(ColumnContent_dom{i})>2
##         for j = 2:length(ColumnContent_dom{i})
##             Sub2Domain_Recover(:,ColumnContent_dom{i}(j)) = lamda_SubOnly(:,i);
##         end
##     end
## end
for i = 1:length(ColumnContent_dom)
  for j = 1:length(ColumnContent_dom{i})
    Sub2Domain_Recover(:, ColumnContent_dom{i}{j}) = lamda_SubOnly(:,i);
  end
end
Sub2Domain_Recover(Sub2Domain_Recover>0.99) = 0.99;

%%-- Predict the Compound-Protein Interactions.
nrow = size(Drug2Sub,1);
ncol = size(Protein2Domain,1)
com2pro = zeros(nrow,ncol);

for i = 1:nrow
  result = cellfun(@GIFTinloop,
                   num2cell(repmat(i,ncol,1)),
                   num2cell(vec2mat(1:ncol,1)),
                   'UniformOutput',false);
    com2pro(i,:) = cell2mat(result);
end
save('tcmdrug2protein_GIFT_songpeng_150410.mat','com2pro');

target_row = importdata('betulinic_target_row.txt'); % selected row
selected_com2pro = com2pro(:,target_row)';
com2pro(:,820)

%% Hierachical Clustering.
%cg1 = clustergram(com2pro,'RowLabels',[tcm2subs(:,1);drug2subs(:,1)],'ColumnLabels',proteinnames,'Standardize','none');



