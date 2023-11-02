function [Mfit,cflag] = postEstCopulaMhlD(parOptBank,Mfit,pIndMap,pOrg)

disp('Forming new copula...')
clear copU xf ilist
% Set Weights from Accuracy
A = (parOptBank(end,:)-1).^-1;


% Compute Mahalanobis Distance
MD = [];
for i = 1:size(parOptBank,2)
    X = parOptBank(pIndMap,i)';
    MD(i) = (X-Mfit.prior.Mu')*inv(Mfit.prior.Sigma)*(X-Mfit.prior.Mu')';
end
MD = 1./MD; % inverse so largest weight the most
W = A+ (MD);
W = W-min(W);
W = W./sum(W);

% First form kernel density estimates for each optimized
% parameter
clear copU
for i = 1:size(pIndMap,1)
    x = parOptBank(pIndMap(i),:); % choose row of parameter values
    bwid(i) = KSDensityCVWidth(x,x,W,[-2 2],25,'cdf');
    copU(i,:) = ksdensity(x,x,'function','cdf','Weights',W,'width',bwid(i)); % KS density estimate per parameter
    xf(i,:) = x;
end
try
    [Rho,nu] = copulafit('t',copU','Method','ApproximateML'); % Fit copula
    % Save outputs that specify the copula
    Mfit.xf = xf;
    Mfit.ks = copU;
    Mfit.nu = nu;
    Mfit.bwid = bwid;
    Mfit.tbr2 = parOptBank(end,1); % get best fit
    Mfit.Pfit = spm_unvec(mean(parOptBank,2),pOrg);
    Mfit.BPfit = spm_unvec(parOptBank(1:end-1,1),pOrg);
    Mfit.Rho = Rho;
    %     Mfit_hist = Mfit;
    %%% Plot posterior, Rho, and example 2D/3D random draws from copulas
%     if R.plot.flag == 1
%         figure(3)
%         clf
%         plotDistChange_KS(Rho,nu,xf,pOrg,pInd,R)
%     end
    cflag = 1;
catch
    disp('Copula estimation is broken, its probably the following:')
    disp('The estimate of Rho has become rank-deficient.  You may have too few data, or strong dependencies among variables.')
    cflag = 0;
end

Mfit.Mu = mean(parOptBank(pIndMap,:),2);
Mfit.Sigma = cov(parOptBank(pIndMap,:)');


if cflag == 1
    
end