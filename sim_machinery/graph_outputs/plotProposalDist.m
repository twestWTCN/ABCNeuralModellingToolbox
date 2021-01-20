function plotProposalDist(R,Mfit,pardraw,pIndMap,cflag)
ls = '-' ;
% Univariate plots
% get priors
preMu  = Mfit.prior.Mu;
preSig = diag(Mfit.prior.Sigma);
cmap = linspecer(5);
xrange = R.SimAn.pOptBound(1):.05:R.SimAn.pOptBound(2);
QL = length(pIndMap);
if QL>5
    QL = 1:ceil(numel(pIndMap)/5):numel(pIndMap); % take the first 5 pars
elseif QL<=5
    QL = 1:numel(pIndMap);
end
% 
subplot(2,2,1)
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
        poln = plotKSPDF(R.SimAn.pOptRange,Mfit.xf(QL(Q),:),r(:,QL(Q)),Mfit.bwid(Q),'-',cmap(Q,:));
    end
end

xlabel('mu')
ylabel('p(mu)')
% ylim([0 1]);
xlim(R.SimAn.pOptBound.*0.2); grid on
title('Approximate Posterior Distribution')
axis square

subplot(2,2,2)
if cflag == 0
    S = Mfit.Sigma; % covariance matrix
    D = sqrt(diag(Mfit.Sigma)); %standard deviations
    Rho = S./(D*D');
else
    Rho = Mfit.Rho;
end
imagesc(Rho);
set(gca,'YDir','normal')
title('Copula Covariance')
set(gca,'XTick',1:size(Rho,1))
set(gca,'YTick',1:size(Rho,1))
set(gcf,'Position',[2.5 617 884 383])
colormap(brewermap(128,'*RdYlBu'))
caxis([-1 1])
axis square
% % Multivariate plots
subplot(2,2,3)
title('2D Draws')
x1 = pardraw(QL(1),:);
y1 = pardraw(QL(2),:);
scatter(x1,y1);
grid on
xlim(R.SimAn.pOptBound.*0.2)
ylim(R.SimAn.pOptBound.*0.2)

xlabel('M1 Time Constant'); ylabel('M1 Synaptic Gain')
set(get(gca,'children'),'marker','.')
axis square

subplot(2,2,4)
title('3D Draws')
z1 = pardraw(QL(3),:);
scatter3(x1,y1,z1)
xlim(R.SimAn.pOptBound.*0.2)
ylim(R.SimAn.pOptBound.*0.2)
zlim(R.SimAn.pOptBound.*0.2)
axis square

xlabel('M1->STR A'); ylabel('STR Gain'); zlabel('M1->STR D');
set(get(gca,'children'),'marker','.')

% set(gcf,'Position',[2.5000  272.0000  903.0000  728.0000])
