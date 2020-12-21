function [r2,pMAP,feat_sim,xsims,xsims_gl,wflag,Rmod] = getMultiStartData(R,intag,multiswitch)
R.out.tag = intag;
R.out.dag = sprintf([R.out.tag '_M%.0f'],1);

[Rmod,m,p,parBank,permMod] = loadABCData_160620(R);
Rmod = setSimTime(Rmod,256);
u = innovate_timeseries(Rmod,m);
pMAP = permMod.MAP;
if multiswitch == 2
    rng(5342)
    x = spm_vec(pMAP);
    [pInd,pMu,pSig] = parOptInds_110817(R,p,m.m); % in structure form
    xMu = x(spm_vec(pMu));
    x(spm_vec(pMu)) = xMu + (sign(randn(size(xMu)))*0.4);
    pMAP= spm_unvec(x,pMAP);
end
[r2,pnew,feat_sim,xsims,xsims_gl,wflag] = computeSimData_160620(Rmod,m,u,pMAP,0,1);
Rmod.data.feat_emp = feat_sim{1};
Rmod.data.feat_xscale{1} = R.frqz;
rng('shuffle')