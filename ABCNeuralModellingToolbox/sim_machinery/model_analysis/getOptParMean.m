function optP = getOptParMean(m,p,R,parBank)
%% Resample parameters
% load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD\outputs\CSD_ABC_neatmodel_iterate2\parBank_CSD_ABC_neatmodel_iterate2_2017818.mat')
% load([R.rootn 'outputs\' R.out.tag '2\parBank_' R.out.tag '2_' d '.mat'])
% figure
% hist(parOptBank(end,:),[-1:.1:1]); xlim([-1 1])
eps = R.analysis.modEvi.eps;
N = R.analysis.modEvi.N;
parOptBank = parBank(1:end-1,parBank(end,:)>eps);
% Ensure manageable size
if size(parOptBank,2)>2^10
    parOptBank = parOptBank(:,1:2^10);
end

% Compute indices of optimised parameter
pInd = parOptInds_110817(R,p,m.m); % in structure form
pIndMap = spm_vec(pInd); % in flat form
R.SimAn.minRank = ceil(size(pIndMap,1)*1.1);
xf = zeros(size(pIndMap,1),size(parOptBank,2));
for i = 1:size(pIndMap,1)
    x = parOptBank(pIndMap(i),:); % choose row of parameter values
    xf(i,:) = x;
end

disp('Drawing from copula...')
r = copularnd('t',R.Mfit.Rho,R.Mfit.nu,N);
clear x1
for Q = 1:size(xf,1)
    x1(Q,:) = ksdensity(xf(Q,:),r(:,Q),'function','icdf');
end
% setup pars from base
clear base
base = repmat(spm_vec(p),1,N);
for i = 1:N
    base(pIndMap,i) = x1(:,i);
end
optP = spm_unvec(mean(base,2),p);