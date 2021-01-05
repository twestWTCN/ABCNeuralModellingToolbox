function [r2,pMAP,feat_sim,xsims,xsims_gl,wflag,Rmod] = getModelData(R,intag,modN)
R.out.tag = intag;
R.out.dag = sprintf([R.out.tag '_M%.0f'],modN);

[Rmod,m,p,parBank,permMod] = loadABCData_160620(R);
Rmod = setSimTime(Rmod,256);
u = innovate_timeseries(Rmod,m);
pMAP = permMod.MAP;
[r2,pnew,feat_sim,xsims,xsims_gl,wflag] = computeSimData_160620(Rmod,m,u,pMAP,0,1);
Rmod.data.feat_emp = feat_sim{1};
Rmod.data.feat_xscale{1} = R.frqz;
