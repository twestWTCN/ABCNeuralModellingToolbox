function [R,m,permMod,xsimMod] = getSimModelData_Draw(Rorg,modID,simtime,allsimchan,nDraws)
% This function will retrieve the parameters from the ABC fit of
% model (modID) and will simulate some data for period specified in simtime for nDraws.

if nargin<4
    allsimchan = 0;
end
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
    R.obs.gainmeth = {'boring'};
    R.obs.trans.norm = 0;
    R.analysis.drawCop = 0;
    R.analysis.drawCoV = 1;
    
    % This allows you to recover all simulated channels, not just those in
    % the empirical data
    if allsimchan == 1
        R.chloc_name = Rorg.chloc_name;
        R.chsim_name = R.chloc_name;
        R.obs.obsstates = Rorg.obs.obsstates;
    end
    %%
    R.analysis.modEvi.N = nDraws; % Number of draws
    R = setSimTime(R,simtime);
    % You need to reduce the parameter space to just 'A', from the original
    % pOptList;
    [R2,pIndMap,pMuMap,pSigMap] = reassignParOpts(R,p,m,{'.A'},5);

    par = postDrawMVN(R2,R2.Mfit,p,pIndMap,pSigMap,R2.analysis.modEvi.N);
    
    % Setup progress Monitor
    if isempty(gcp('nocreate'))
        parpool
    end
    ppm = ParforProgMon('Draw in progress', nDraws, 1, 400, 100);
    
    % Simulate Models
    parfor jj = 1:numel(par)
        pnew = par{jj};
        u = innovate_timeseries(R,m);
        u{1} = u{1}.*sqrt(R.IntP.dt);
        [r2,pnew,feat_sim,xsims,xsims_gl,wflag] = computeSimData(R,m,u,pnew,0);
        
        %     R.plot.outFeatFx({},{feat_sim},R.data.feat_xscale,R,1)
        wfstr(jj) = any(wflag);
        r2rep{jj} = r2;
        par_rep{jj} = pnew;
        feat_rep{jj} = feat_sim;
        xsims_rep{jj} = xsims_gl;
        disp(jj)
        ppm.increment();
    end
    permMod.wflag= wfstr;
    permMod.r2rep = r2rep;
    permMod.par_rep = par_rep;
    permMod.feat_rep = feat_rep;
    xsimMod = xsims_rep;
    
end