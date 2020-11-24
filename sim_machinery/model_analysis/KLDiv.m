function [KL,DKL,R,Mfit] = KLDiv(R,Mfit,p,m,npflag)
%% Resample parameters
% Compute indices of optimised parameter
pVec = spm_vec(p);
pInd = parOptInds_110817(R,p,m.m); % in structure form
pIndMap = spm_vec(pInd); % in flat form

precInd = parOptIndsPrec_110817(R,p,m.m,2);
precIndMap = spm_vec(precInd);
p = pVec(pIndMap);

% R.SimAn.minRank = ceil(size(pIndMap,1)*1.1);
% xf = zeros(size(pIndMap,1),size(parOptBank,2));
% for i = 1:size(pIndMap,1)
%     xf(i,:) = parOptBank(pIndMap(i),:); % choose row of parameter values
% end

priorMu = Mfit.prior.Mu;
priorSig = Mfit.prior.Sigma;

incpar = find(diag(priorSig)>1e-8);

mu1 = priorMu(incpar);
sig1 = priorSig(incpar,incpar);

KL = [];
if npflag == 1
    xf = Mfit.xf;
    r = copularnd('t',Mfit.Rho,Mfit.nu,R.analysis.modEvi.N);
    for Q = incpar'
        %     mu_prior(Q) = mean(r(:,Q));
        %     mu_post(Q) = priorVec(pIndMap(Q));
        
        y = ksdensity(xf(Q,:),r(:,Q),'function','icdf');
        [y,f] = ksdensity(y,R.SimAn.pOptRange);
        %     x1
        %     plot(f,x1,'r')
        %     hold on
        x1 = normpdf(f,priorMu(Q),priorSig(Q,Q));
        %     plot(f,y,'b')
        
        % Add alpha to get rid of zeros
        y = y + 1e-10;
        x1 = x1 + 1e-10;
        KL(Q) = -sum(y.*log(x1./y));
        %     plot(f,x1./y,'g')
    end
    
    
    sig2 = cov(xf(incpar,:)');
    mu2 = mean(xf(incpar,:)')';
    
else
    sig2 = Mfit.Sigma;
    mu2 = Mfit.Mu;
end


DKL = 0.5*(trace(sig2\sig1) + ((mu2-mu1)'*sig2)'\(mu2-mu1) - size(mu1,1) + log(det(sig2)/det(sig1)));
Mfit.DKL = DKL;
Mfit.KL = KL;
Mfit.npflag = npflag;
R.Mfit = Mfit;
