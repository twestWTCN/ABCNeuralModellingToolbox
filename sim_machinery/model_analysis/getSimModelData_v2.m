function [R,m,permMod,xsimMod] = getSimModelData_v2(R,modID,simtime)
% This function will specficailly get the parameters from postieror fit of
% model and will simulate some data.
for pn = 1:numel(modID)
    rng(2223) % Ensure random elements are the same
    % Load Model Data
    R.out.dag = sprintf(['NPD_' R.out.tag '_M%.0f'],modID(pn));
    [R,m,p,parBank] = loadABCData(R);
    a = eval(['@MS_rat_' R.out.tag '_Model' num2str(modID(pn))]);
    [dum prior] = a(R);
    R.Mfit.prior = prior;
    R.analysis.modEvi.eps = parBank(end,R.SimAn.minRank);
    parOptBank = parBank(1:end,1:2^10);
    R.parOptBank = parOptBank;
    R.obs.gainmeth = {};
    R.obs.trans.norm = 0;
    R.analysis.modEvi.N = 2000;
    R.analysis.BAA.flag = 1;
    R.analysis.BAA.redmeth = 'UQ'; % average samples to get parameters
    R = setSimTime(R,simtime);
       
    % With Hyperdirect
    [permMod{pn} xsimMod{pn}] = modelProbs(m.x,m,p,R);
end