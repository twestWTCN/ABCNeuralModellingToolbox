function [R,m,permMod,xsimMod] = getSimModelData_v3(Rorg,modID,simtime,allsimchan)
if nargin<4
    allsimchan = 0;
end
% This function will  get the parameters from the ABC fit of
% model (modID) and will simulate some data for period simtime.
for pn = 1:numel(modID)
    rng(2223) % Ensure random elements are the same
    % Load Model Data
    Rorg.out.dag = sprintf(['NPD_' Rorg.out.tag '_M%.0f'],modID(pn));
    [R,m,p,parBank] = loadABCData(Rorg);
    a = eval(['@MS_rat_' Rorg.out.tag '_Model' num2str(modID(pn))]);
    [dum prior] = a(R);
    R.Mfit.prior = prior;
    R.analysis.modEvi.eps = parBank(end,R.SimAn.minRank);
    parOptBank = parBank(1:end,1:2^10);
    R.parOptBank = parOptBank;
    R.obs.gainmeth = {};
    
    % This allows you to recover all simulated channels, not just those in
    % the empirical data
    if allsimchan == 1
        R.chloc_name = Rorg.chloc_name;
        R.chsim_name = R.chloc_name;
        R.obs.obsstates = Rorg.obs.obsstates;
    end
    %%
    R.analysis.modEvi.N = 2000; % Number of draws
    R.analysis.BAA.flag = 1;
    R.analysis.BAA.redmeth = 'UQ'; % average samples to get parameters
    R = setSimTime(R,simtime);
    
    % Legacy Models Require Mfit fields in the .prior field (copy over)
    R.Mfit.prior = R.Mfit;
    warning('getSimModelData_v3 is using legacy setup!')
    % With Hyperdirect
    [permMod{pn} xsimMod{pn}] = modelProbs(m.x,m,p,R);
end