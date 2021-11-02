function [R,m,permMod,xsimMod] = getSimModelData_160620(Rorg,modID,simtime,allsimchan)
if nargin<4
    allsimchan = 0;
end
% This function will  get the parameters from the ABC fit of
% model (modID) and will simulate some data for period simtime.
for pn = 1:numel(modID)
    rng(2223) % Ensure random elements are the same
    % Load Model Data
    Rorg.out.dag = sprintf([Rorg.out.tag '_M%.0f'],modID);
    [R,m,p,parBank] = loadABCData_160620(Rorg);
    R.analysis.modEvi.eps = parBank(end,R.SimAn.minRank);
    parOptBank = parBank(1:end,1:2^10);
    R.parOptBank = parOptBank;
    R.obs.gainmeth = {};
    
    R.analysis.modEvi.N = 1; % Number of draws
    R.analysis.BAA.flag = 1;
    R.analysis.BAA.redmeth = 'average'; % average samples to get parameters
    R = setSimTime(R,simtime);
    
    % Legacy Models Require Mfit fields in the .prior field (copy over)
    R.Mfit.prior = R.Mfit;
    warning('getSimModelData_v3 is using legacy setup!')
    % With Hyperdirect
    [permMod{pn} xsimMod{pn}] = modelProbs_160620(R,m.x,m,p,R);
end

R.path = Rorg.path;