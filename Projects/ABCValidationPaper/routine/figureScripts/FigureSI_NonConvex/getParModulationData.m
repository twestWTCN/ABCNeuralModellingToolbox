

function [r2,pMAP,feat_sim,xsims,xsims_gl,wflag,Rmod] = getParModulationData(R,intag,field2mod)
R.out.tag = intag;
R.out.dag = sprintf([R.out.tag '_M%.0f'],10);

[Rmod,m,p,parBank,permMod] = loadABCData_160620(R);
Rmod = setSimTime(Rmod,32);
u = innovate_timeseries(Rmod,m);
pMAP = permMod.MAP;

A = linspace(-0.25,0.25,25);
B = linspace(-1,1,25);
r2 = zeros(numel(A),numel(B));
for i = 1:numel(A)
    r2tmp = ones(1,numel(B));
    ptmp_i = pMAP;
    ptmp_i.D(4,1) = A(i);
    parfor j = 1:numel(B)
        ptmp_ij = ptmp_i;
        ptmp_ij.A{1}(1,6) = B(j);
        r2tmp(j) = computeSimData_160620(Rmod,m,u,ptmp_ij,0,1);
        disp([i,j])
    end
    r2(i,:) = r2tmp;
end

a = [];