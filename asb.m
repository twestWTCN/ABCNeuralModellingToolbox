

[~,m] = loadABCData_160620(R);

[pInd,pMu,pSig] = parOptInds_110817(R,MAP(1),2); % in structure form
pMuMap = spm_vec(pMu);
pSigMap = spm_vec(pSig);
parNamesFull = getParFieldNames(MAP(1),m);
parNameMu = parNamesFull(pMuMap);
parNameSig = parNamesFull(pSigMap);


Xbar = []; sig = [];
for M = 1:10
a = spm_vec(MAP(M));
Xbar(:,M) = a(spm_vec(pMu));
sig(:,M) = a(spm_vec(pSig));
end
b = bar(Xbar)
a = gca;
a.XTickLabel = parNameMu
hold on
errorbar(1:size(Xbar,1),mean(Xbar,2),mean(sig,2),'LineStyle','none')