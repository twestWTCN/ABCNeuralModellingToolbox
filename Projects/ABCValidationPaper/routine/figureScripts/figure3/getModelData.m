function [r2,pMAP,feat_sim,xsims,xsims_gl,wflag,Rmod] = getModelData(R,intag,modN,fixNoise,dt)
if nargin<4
    fixNoise =0;
end
R.out.tag = intag;
R.out.dag = sprintf([R.out.tag '_M%.0f'],modN);

[Rmod,m,p,parBank,permMod] = loadABCData_160620(R);
if nargin == 5
    Rmod.IntP.dt = dt;
    Rmod.obs.trans.norm 
end
Rmod = setSimTime(Rmod,256);
if fixNoise
    rng(6433)
end
u = innovate_timeseries(Rmod,m);
pMAP = permMod.MAP;
[r2,pnew,feat_sim,xsims,xsims_gl,wflag] = computeSimData_160620(Rmod,m,u,pMAP,0,1);
Rmod.data.feat_emp = feat_sim{1};
% Rmod.data.feat_xscale{1} = R.frqz;
