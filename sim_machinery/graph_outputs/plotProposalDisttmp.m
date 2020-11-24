function plotProposalDist(R,Mfit,pIndMap,cflag)
ls = '-' ;
% Univariate plots
% get priors
preMu  = Mfit.prior.Mu;
preSig = diag(Mfit.prior.Sigma);
cmap = linspecer(5);
xrange = R.SimAn.pOptBound(1):.05:R.SimAn.pOptBound(2);
QL = length(pIndMap);
if QL>3
    QL = 1:ceil(numel(pIndMap)/3):numel(pIndMap); % take the first 5 pars
end

if cflag == 0 % Normal Dist
    % get proposal
    proMu  = Mfit.Mu;
    proSig = diag(Mfit.Sigma);
    
    for Q = 1:numel(QL)
        peln = plotNormalPDF(xrange,preMu(QL(Q)),preSig(QL(Q)),'--',cmap(Q,:));
        poln = plotNormalPDF(xrange,proMu(QL(Q)),proSig(QL(Q)),'-',cmap(Q,:));
    end
elseif cflag == 1 % Kernel Density Estimates
    ls = '-';
    r = copularnd('t',Mfit.Rho,Mfit.nu,500);
    for Q = 1:numel(QL)
        peln = plotNormalPDF(xrange,preMu(QL(Q)),preSig(QL(Q)),'--',cmap(Q,:));
        poln = plotKSPDF(R.SimAn.pOptRange,Mfit.xf(QL(Q),:),r(:,QL(Q)),'-',cmap(Q,:));
    end
end

xlabel('\mu')
ylabel('p(\mu)')
% ylim([0 1]);
xlim(R.SimAn.pOptBound.*0.2)
title('Approximate Posterior Distribution')

% subplot(2,2,2)
% imagesc(Rho)
% set(gca,'YDir','normal')
% title('Copula Covariance')
% set(gca,'XTick',1:size(Rho,1))
% set(gca,'YTick',1:size(Rho,1))
% set(gcf,'Position',[2.5 617 884 383])
% % Multivariate plots
% 
% r = copularnd('t',Rho,nu,1000);
% subplot(2,2,3)
% title('2D Sample Drawn from Copula')
% if strcmp(R.projectn,'MVAR')
%     i = pInd.params(1);
%     j =pInd.params(2);
% else
%     i = spm_vec(pInd);
%     i = i(4);
%     j = spm_vec(pInd);
%     j = j(5);
% end
% i = find(indFlat==i);
% j = find(indFlat==j);
% u1 = r(:,i);
% v1 = r(:,j);
% x1 = ksdensity(xf(i,:),u1,'function','icdf');
% y1 = ksdensity(xf(j,:),v1,'function','icdf');
% scatter(x1,y1);
% xlim(R.SimAn.pOptBound.*0.2)
% ylim(R.SimAn.pOptBound.*0.2)
% 
% xlabel('M1 Time Constant'); ylabel('M1 Synaptic Gain')
% set(get(gca,'children'),'marker','.')
% 
% subplot(2,2,4)
% title(' 3D Sample Drawn from Copula')
% if strcmp(R.projectn,'MVAR')
%     i = pInd.params(1);
%     j =pInd.params(2);
%     k = pInd.params(3);
% else
%     i = spm_vec(pInd);
%     i = i(4);
%     j = spm_vec(pInd);
%     j = j(5);
%     k = spm_vec(pInd);
%     k = k(2);
% end
% i = find(indFlat==i);
% j = find(indFlat==j);
% k = find(indFlat==k);
% 
% u1 = r(:,i);
% v1 = r(:,j);
% w1 = r(:,k);
% x1 = ksdensity(xf(i,:),u1,'function','icdf');
% y1 = ksdensity(xf(j,:),v1,'function','icdf');
% z1 = ksdensity(xf(k,:),w1,'function','icdf');
% 
% scatter3(x1,y1,z1)
% xlim(R.SimAn.pOptBound.*0.2)
% ylim(R.SimAn.pOptBound.*0.2)
% zlim(R.SimAn.pOptBound.*0.2)
% 
% xlabel('M1->STR A'); ylabel('STR Gain'); zlabel('M1->STR D');
% set(get(gca,'children'),'marker','.')
% 
% set(gcf,'Position',[2.5000  272.0000  903.0000  728.0000])
