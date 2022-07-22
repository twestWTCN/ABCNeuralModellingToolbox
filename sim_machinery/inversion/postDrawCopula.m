function [par,MAP] = postDrawCopula(R,Mfit,pOrg,pIndMap,pSigMap,rep,permScale)
if nargin<7
    permScale = 0;
end
disp('Drawing from copula...')
xrange = R.SimAn.pOptBound(1):.005:R.SimAn.pOptBound(2);

r = copularnd('t',Mfit.Rho,Mfit.nu,rep);
clear x1
xf = Mfit.xf;
for Q = 1:size(xf,1)
    %     bwid(Q) = KSDensityCVWidth(xf(Q,:),r(:,Q),repmat(1./size(xf,1),1,size(xf,2)),[-2 3],25,'icdf');
    x1(Q,:) = ksdensity(xf(Q,:),r(:,Q),'function','icdf','width',Mfit.bwid(Q));
    % permute samples
    x1(Q,:) = x1(Q,:) + (permScale*std(x1(Q,:))*randn(size(x1(Q,:))));
    
    % Construct MAP estimates from peak of marginal distribution
    [px,fx] = ksdensity( x1(Q,:),xrange,'function','pdf','width',Mfit.bwid(Q));
    mapq(Q) = fx(px==max(px));
end
% setup pars from base
clear base
base = repmat(spm_vec(pOrg),1,rep);
for i = 1:rep
    base(pIndMap,i) = x1(:,i);
    base(pSigMap,i) = diag(Mfit.Sigma); % variance is the diagonal of the covariance matrix
    par{i} = spm_unvec(base(:,i),pOrg);
end
baseSave = base;
% Set up MAP estimate
base = spm_vec(pOrg);
MAP = base;

MAP(pIndMap) = mapq;             % this is the mode of the density
MAP(pSigMap) = diag(Mfit.Sigma); % this is the variance of the estimate (derived from samples, not the density)
MAP = spm_unvec(MAP,pOrg);

if R.plot.flag == 1
    if ~ishandle(3)
        figure(3)
    else
        set(groot,'CurrentFigure',3);  clf
    end
    clf
    cflag= 1;
    pardraw = baseSave(pIndMap,:);
    plotProposalDist(R,Mfit,pardraw,pIndMap,cflag)
end

