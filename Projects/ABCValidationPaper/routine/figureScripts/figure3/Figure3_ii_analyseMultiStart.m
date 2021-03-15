function Figure3_ii_analyseMultiStart(R)
%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3- (II) Multistart Analysis
%%%%%%%%%%%%%%%%%%%%%%%%
R.out.tag = 'figure3_MultiStart'; % Task tag
R = ABCsetup_partI_STNGPe(R);

modelspec = eval(['@MS_rat_STN_GPe_ModComp_Model' num2str(1)]);
[R,p,m] = modelspec(R);
N = 8; % Number of multistarts

% Load Priors
load([R.path.rootn '\outputs\' R.path.projectn '\figure3_MultiStart\MultiStartDataFeatures.mat'])


[pInd,parMu,parSigMap] = parOptInds_110817(R,p,m.m); % in structure form
% Set Par Names
parNames = getParFieldNames(p,m);

parSel = 1:12;
% Form parameter indicies
pMuMap = spm_vec(parMu); % in flat form
pMuMap = pMuMap(parSel);
pSigMap = spm_vec(parSigMap); % in flat form
pSigMap = pSigMap(parSel);
parNames = parNames(pMuMap);

% Get priors
for ds = 1:2
    pAct = spm_vec(pMAP{ds});
    pMuAct{ds} = pAct(pMuMap);
    pSigAct{ds} = pAct(pSigMap);
%     A = 1-( pSigAct{ds}.^(1/2));
    pMuAct{ds} = pMuAct{ds}; %.*A;
    pMuAct{ds}(end-1:end) = 0;
end

% Run through multistarts
Inds(1,1) = 0; Inds(2,1) = 0;
for multiStart = 1:(2*N)
    R.out.dag = sprintf('NPD_STN_GPe_MultiStart_M%.0f',multiStart); % 'All Cross'
    [Rout,m,p,parBank,~,parHist,bankSave,kldHist] = loadABCData_160620(R);
    
    
    % parameter history
    parTT = []; r2 = [];
    for i = 1:size(parHist,2)
        parTT(:,i) = spm_vec(parHist(i));
        r2(i) = mean(bankSave{i});
    end
    parT = parTT(pMuMap,:);
    parSig{multiStart} = parTT(pSigMap,:);
    
    
    % Precision
    A = 1-( parSig{multiStart}.^(1/2));
    
    % Save variables
    wParT = parT.*A;
    parWeighted{multiStart} = wParT;
    parMS{multiStart} = parT;
    parConv(:,multiStart) = wParT(:,end);
    parSigConv(:,multiStart) = parTT(pSigMap,end);
    Inds(:,multiStart+1) = [Inds(2,multiStart)+1, Inds(2,multiStart) + size(parT,2)];
    R2ms(multiStart) = r2(end);
    R2track{multiStart} = r2;
    Its(multiStart) = size(parT,2);
    kldTrack{multiStart} = kldHist;
    if multiStart<=N
        pwpe{multiStart} = (wParT-pMuAct{1}).^2;
    elseif multiStart>N
        pwpe{multiStart} = (wParT-pMuAct{1}).^2;
    end
end
Inds(:,1) = [];

mkdir([R.path.rootn '\outputs\' R.path.projectn '\figure3_MultiStart'])
save([R.path.rootn '\outputs\' R.path.projectn '\figure3_MultiStart\MSAsave1.mat'])
% save('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD\routine\rat_STN_GPe\Model_Validation\MultiStartAnalysis\MSAsave3.mat')