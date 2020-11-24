close all; clear
R = simannealsetup_NPD_STN_GPe;
R.out.dag = sprintf('NPD_STN_GPe_MultiStart_M%.0f',1); % 'All Cross'

modelspec = eval(['@MS_rat_STN_GPe_ModComp_Model' num2str(1)]);
[R,p,m] = modelspec(R);
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' R.out.dag '.mat'])
m = varo;

[pInd,parMu,parSigMap] = parOptInds_110817(R,p,m.m); % in structure form
% [pInd] = createParNameField(R,p,m.m); % in structure form
% Set Par Names
parNames = {'GPe \tau';
    'GPe \gamma';
    'STN \tau';
    'STN \gamma';
    'GPe C';
    'STN C';
    'STN \rightarrow GPe';
    'GPe \rightarrow STN';
    'STN D\rightarrow GPe';
    'GPe D\rightarrow STN'; 
    };
parSel = [1 3 5:10];
% Form descriptives
pIndMap = spm_vec(parMu); % in flat form
pIndMap = pIndMap(parSel);
pSigMap = spm_vec(parSigMap); % in flat form
pSigMap = pSigMap(parSel);
parNames = parNames(parSel);
Inds(1,1) = 0; Inds(2,1) = 0;
for multiStart = 1:20
    R.out.dag = sprintf('NPD_STN_GPe_MultiStart_M%.0f',multiStart); % 'All Cross'
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
    Mfit = varo;
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parHist_' R.out.tag '_' R.out.dag '.mat'])
    parHist = varo;
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\bankSave_' R.out.tag '_' R.out.dag '.mat'])
    bankSave = varo;
    
    parTT = []; r2 = [];
    for i = 1:size(parHist,2)
        parTT(:,i) = spm_vec(parHist(i));
        r2(i) = mean(bankSave{i});
    end
    parT = parTT(pIndMap,:);
    parSig{multiStart} = parTT(pSigMap,:);
    %     parSig = parT(pIndMap,:);
    %     A = 1./parSig;
    %     A = A./sum(A,1);
    A = 1-( parSig{multiStart}.^(1/2));
    %     A = 1/(parSig{multiStart}./parSig{multiStart}(:,1));
    %     A(A<0) = 0
    %     A = 1;
    wParT = A.*parT;
    parWeighted{multiStart} = wParT;
    parMS{multiStart} = parT;
    parConv(:,multiStart) = wParT(:,end);
    Inds(:,multiStart+1) = [Inds(2,multiStart)+1, Inds(2,multiStart) + size(parT,2)];
    R2ms(multiStart) = r2(end);
    R2track{multiStart} = r2;
    Its(multiStart) = size(parT,2);
end
Inds(:,1) = [];
T = [parWeighted{:}];

% T1 = [parWeighted{1:10}];
% T2 = [parWeighted{11:20}];
% 
% [A,B,r,U,V,stats] = canoncorr(T1,T2) 

D = pdist(T','euclidean');
[Y,eigvals] = cmdscale(squareform(D));

save('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD\routine\rat_STN_GPe\Model_Validation\MultiStartAnalysis\MSAsave.mat')