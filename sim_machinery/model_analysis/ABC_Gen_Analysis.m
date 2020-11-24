function [] = ABC_Gen_Analysis(R,p,m,parBank)

% Create parameter index map
pInd = parOptInds_110817(R,p,m.m,2); % in structure form
pIndMap = spm_vec(pInd); % in flat form


% Re-compute optimised parbank
parOptBank = parBank(:,parBank(end,:)>R.analysis.modEvi.eps);
R.parOptBank = parOptBank;

for i = 1:size(pIndMap,1)
    x = R.parOptBank(pIndMap(i),:); % choose row of parameter values
    copU(i,:) = ksdensity(x,x,'function','cdf'); % KS density estimate per parameter
    xf(i,:) = x;
end

[Rho,nu] = copulafit('t',copU','Method','ApproximateML'); % Fit copula
% Save outputs that specify the copula
R.Mfit.Rho = Rho;
R.Mfit.xf = xf;
R.Mfit.nu = nu;


%% Triangle plot of MV posterior distribution
% figure(1)
% triangle_plot_mdist(R,p,m,xf)

%% Compute Model Evidence
gcp
figure(2);
% R.analysis.modEvi.N = 128;
modelProbs(m.x,m,p,R)

%% Plot Simulated Features with Confidence intervals
figure(3)
PlotFeatureConfInt_gen(R)
saveallfiguresFIL_n([R.rootn '\analysis\' R.out.tag '1\featConfInts'],'-jpg',1,'-r200',2);

optP = getOptParMean(m,p,R,d);

%% Parameter Sweep across STN/GPe Subcircuit Connections
parsweep.R = '.A{1}(3,4)';
parsweep.Rname = 'STN-> GPe Connections Strength';
parsweep.Q = '.A{2}(4,3)';
parsweep.Qname = 'GPe-|STN Connections Strength';

% Conduct sweep
parsweep = modelBetaParSweep(m,optP,parsweep,R);
pathstr = [R.rootn 'analysis\parsweeps\'];
mkdir(pathstr);
save([pathstr '\A_STNGPe_AGPeSTN_parsweep'],'parsweep');

% Plot heatmap
figure(4)
parSweepPlot(R,parsweep)
 saveallfiguresFIL_n([R.rootn '\analysis\' R.out.tag '\parSweepSTNGPe.jpg'],'-jpg',1,'-r200',2);
