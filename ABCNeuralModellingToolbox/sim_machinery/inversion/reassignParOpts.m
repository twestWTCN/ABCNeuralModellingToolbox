function [Rout,pIndMap,pMuMap,pSigMap] = reassignParOpts(Rin,p,m,pOptList,rscale)
% This function will crop the existing parameter structure to the reduced
% list given in pOptList
    [pInd,pMu,pSig] = parOptInds_110817(Rin,p,m.m); % in structure form
    pIndlistFull = spm_vec(pInd);
    
    % Now find just for A
    Rtmp.SimAn.pOptList = pOptList;
    [pInd,pMu,pSig] = parOptInds_110817(Rtmp,p,m.m); % in structure form
    pIndlistA = spm_vec(pInd);
    
    for i = 1:numel(pIndlistA)
        Aind(i) = find(pIndlistFull==pIndlistA(i));
    end
    
    % Form descriptives
    pIndMap = spm_vec(pInd); % in flat form
    pMuMap = spm_vec(pMu);
    pSigMap = spm_vec(pSig);
    Rin.plot.flag = 0;
    
    Rout = Rin;
    Rout.Mfit.Sigma = Rout.Mfit.Sigma(Aind,Aind);
    I = sub2ind(size(Rout.Mfit.Sigma),1:numel(Aind),1:numel(Aind));
    Rout.Mfit.Sigma(I) = Rout.Mfit.Sigma(I)*rscale;
    Rout.Mfit.Mu = Rout.Mfit.Mu(Aind);