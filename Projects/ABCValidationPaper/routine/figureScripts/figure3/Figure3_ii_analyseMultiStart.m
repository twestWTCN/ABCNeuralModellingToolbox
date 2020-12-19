close all; clear
%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3- (II) Multistart Analysis
%%%%%%%%%%%%%%%%%%%%%%%%
%This should link to your repo folder
repopath = 'C:\Users\timot\Documents\GitHub\ABCNeuralModellingToolbox';
%This should be your projectname
projname = 'ABCValidationPaper';
R = ABCAddPaths(repopath,projname);


R.out.tag = 'figure3_MultiStart'; % Task tag
R = ABCsetup_partI_STNGPe(R);

modelspec = eval(['@MS_rat_STN_GPe_ModComp_Model' num2str(1)]);
[R,p,m] = modelspec(R);
N = 10; % Number of multistarts


[pInd,parMu,parSigMap] = parOptInds_110817(R,p,m.m); % in structure form
% Set Par Names
parNames = getParFieldNames(p,m);

parSel = 1:12;
% Form descriptives
pMuMap = spm_vec(parMu); % in flat form
pMuMap = pMuMap(parSel);
pSigMap = spm_vec(parSigMap); % in flat form
pSigMap = pSigMap(parSel);
parNames = parNames(pMuMap);
Inds(1,1) = 0; Inds(2,1) = 0;
for multiStart = 1:(2*N)
    R.out.dag = sprintf('NPD_STN_GPe_MultiStart_M%.0f',multiStart); % 'All Cross'
    [Rout,m,p,parBank,~,parHist,bankSave,kldHist] = loadABCData_160620(R);

    
    parTT = []; r2 = [];
    for i = 1:size(parHist,2)
        parTT(:,i) = spm_vec(parHist(i));
        r2(i) = mean(bankSave{i});
    end
    parT = parTT(pMuMap,:);
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
    kldTrack{multiStart} = kldHist;
end
Inds(:,1) = [];
T = [parWeighted{:}];

% T1 = [parWeighted{1:10}];
% T2 = [parWeighted{11:20}];
% 
% [A,B,r,U,V,stats] = canoncorr(T1,T2) 


D = pdist(T','euclidean');
[Y,eigvals] = cmdscale(squareform(D));
save([R.path.rootn '\outputs\' R.path.projectn '\MultiStartAnalysis\MSAsave1.mat'])
% save('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD\routine\rat_STN_GPe\Model_Validation\MultiStartAnalysis\MSAsave3.mat')