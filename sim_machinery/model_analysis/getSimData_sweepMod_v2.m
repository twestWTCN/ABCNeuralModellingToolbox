function [R,permMod,xsimMod] = getSimData_sweepMod_v2(R,modID,simtime)

for pn = 1:numel(modID)
    % Load Model
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' R.out.dag '.mat'])
    m = varo;
    % load modelfit
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
    A = varo;
    p = A.BPfit;
    % Load Options
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\R_' R.out.tag '_' R.out.dag '.mat'])
    R = varo;
    
    R.Mfit = A;
    a = eval(['@MS_rat_' R.out.tag '_Model' num2str(modID(pn))]);
    [dum prior] = a(R);
    R.Mfit.prior = prior;
    % load parbank?
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parBank_' R.out.tag '_' R.out.dag '.mat'])
    parBank =  varo;
    
    R.analysis.modEvi.eps = parBank(end,R.SimAn.minRank);
    parOptBank = parBank(1:end-1,2^10);
    
    R.parOptBank = parOptBank;
    R.obs.gainmeth = R.obs.gainmeth(1);
    
    R.analysis.modEvi.N = 2000;
    R.analysis.BAA = 1;
    R = setSimTime(R,simtime);
    
    % With Hyperdirect
    [permMod{pn} xsimMod{pn}] = modelProbs(m.x,m,p,R);
end