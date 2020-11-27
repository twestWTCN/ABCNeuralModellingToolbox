function [R,parBank] = SimAn_ABC_201120(R,p,m,parBank)
%%%% APROXIMATE BAYESIAN COMPUTATION for
%%%% HIGH DIMENSIONAL DYNAMICAL MODELS
% ---- 25/03/20---------------------------
% This annealing function uses approximate Bayesian computation in order to
% estimate the posterior parameter distributions, using a shifting epsilon
% that moves with a cooling schedule.
% Notes
% This version collapses parameter resamples into a single function
% adapted for generic usage
% ic - starting conditions
% u - external input (can be empty
% p - structure of parameters
% m - structure of model specfications
% R - settings for annealing, plotting, integration etc.
%
% This function will output R at each annealing loop (assigned into 'base'
% workspace. The field R.mfit is appended with the Rho (correlation matrix) and nu
% (degrees of freedom) of the estimated copula. Multivariate Gaussian
% is also estimated and specified in R.mfit.Mu and R.mfit.Sigma.
% There are several plotting functions which will track the progress of the
% annealing.
% TO DO:
%
%
%
% Timothy West (2018) - UCL CoMPLEX
% / UCL, Wellcome Trust Centre for Human Neuroscience
%%%%%%%%%%%%%%%%%%%%%%
warning('off', 'MATLAB:MKDIR:DirectoryExists');
%% Setup for annealing
if nargin<4
    parBank = [];
end
if ~isfield(R.plot,'flag')
    R.plot.flag = 1; % plotting is default behaviour
end
if isempty(m)
    m.m = 1;
end
if ~isfield(R.plot,'updateflag')
    R.plot.updateflag = 0;
end
pOrg = p; % Record prior parameters.

% Set Fixed Annealing Parameters
eps_prior = -200; % prior eps (needed for gradient approximation);
eps_exp = -12;
eps_act = eps_prior;
delta_act = 0.05;
% Compute indices of parameters to be optimized
[pInd,pMu,pSig] = parOptInds_110817(R,p,m.m); % in structure form
% Form descriptives
pIndMap = spm_vec(pInd); % in flat form
pMuMap = spm_vec(pMu);
pSigMap = spm_vec(pSig);
R.SimAn.minRank = ceil(size(pIndMap,1)*2); %Ensure rank of sample is large enough to compute copula

% set initial batch of parameters from gaussian priors
if isfield(R,'Mfit')
    rep =  R.SimAn.rep(1);
    par = postDrawCopula(R.Mfit,p,pIndMap,rep);
else
    rep = R.SimAn.rep(1);
    ptmp = spm_vec(p);
    Mfit.Mu = ptmp(pMuMap);
    Mfit.Sigma = diag(ptmp(pSigMap).*R.SimAn.jitter);
    Mfit.prior = Mfit;
    par = postDrawMVN(R,Mfit,pOrg,pIndMap,pSigMap,rep);
end
R.Mfit.DKL = 0; % Divergence is zero to begin
parPrec(:,1) = diag(Mfit.Sigma);
itry = 0; cflag = 0;
ii = 1; parOptBank = [];
%% Main Annealing Loop
while ii <= R.SimAn.searchMax
    %% Batch Loop for Replicates for Generation of Pseudodata
    % This is where the heavy work is done. This is run inside parfor. Any
    % optimization here is prime.
    clear xsims_rep feat_sim_rep featbank ACCbank
    ji = 0;
    parnum = (6*2);
    samppar = {}; ACCbank = []; featbank = [];
    while ji < floor(rep/parnum)
       parfor jj = 1:parnum % Replicates for each temperature
            % Get sample Parameters
            parl = (ji*parnum) + jj;
            pnew = par{parl};
            %% Simulate New Data
            u = innovate_timeseries(R,m);
            [r2,pnew,feat_sim] = computeSimData_160620(R,m,u,pnew,0,0);
            % Adjust the score to account for set complexity
            
            [ACC,R2w] = computeObjective(R,r2);
            r2rep(jj) = R2w;
            ACCrep(jj) = ACC;
            par_rep{jj} = pnew;
            %         xsims_rep{jj} = xsims_gl; % This takes too much memory: !Modified to store last second only!
            feat_sim_rep{jj} = feat_sim;
            
%             fprintf(1,'\b\b%.0f',jj/parnum);
        end % End of batch replicates
        if rem(ji,4)
        disp(['Batch ' num2str(ji) ' proposal ' num2str(ii)])
        end
        % Retrieve fits
        
        r2loop = r2rep; %ACCrep;
        % Delete failed simulations
        r2loop(r2loop==1) = -inf;
        r2loop(isnan(r2loop)==1) = -inf;
        r2loop(imag(r2loop)==1) = -inf;
        r2loop(isinf(r2loop)==1) = -inf;
        % Append succesful replicates to bank of params and fits
        %(parameter table, with fits)
        for i = 1:numel(r2loop)
            if ~isinf(r2loop(i))
                parI(:,i) = [full(spm_vec(par_rep{i})); r2loop(i)]';
                parBank = [parBank parI(:,i) ];
            end
        end
        
        % Save data features (for plotting reasons only)
        featbank{ji+1} = feat_sim_rep;
        ACCbank(:,ji+1) = ACCrep;
        samppar{ji+1} = par_rep;
        % Clip parBank to the best (keeps size manageable
        if ~isempty(parBank)
            [dum V] = sort(parBank(end,:),'descend');
            if size(parBank,2)>2^13
                parBank = parBank(:,V(1:2^12));
            else
                parBank = parBank(:,V);
            end
            ACClocbank = computeObjective(R,parBank(end,:));
            parOptBank = parBank(:,ACClocbank>eps_exp);
            if size(parOptBank,2)> R.SimAn.minRank-1
                break
            else
                ji = ji+1;
            end
        else
            ji = ji+1;
        end
        
    end
    
    %% Find the best draws (for plotting only)
    bestfeat = [];
    [b i] = maxk(ACCbank(:),12);
    [jj_best, ji_best] = ind2sub(size(ACCbank),i);
    
    % Simulate best data (plotting outside of parfor)
    pnew = samppar{ji_best(1)}{jj_best(1)};
    u = innovate_timeseries(R,m);
    [~,~,~,~,xsims_gl_best] = computeSimData_160620(R,m,u,pnew,0,0);
    
    for L = 1:numel(i)
        bestfeat{L}{1} = featbank{ji_best(L)}{jj_best(L)}{1};
    end
    bestr2(ii) =  ACCbank(jj_best(1),ji_best(1));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PARAMETER OPTIMIZATION BEGINS HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % icop(1)>itry: Try to form copula with Annealing or Percentile eps
    % icop(2)>itry>icop(1): Find eps to form minimum rank from parbank
    % icop(2)<itry: Try to force
    
    %% Concatanate Batch Results and Decide Acceptance Level Epsilon
    
    %% Find error threshold for temperature (epsilon) and do rejection sampling
    A = parOptBank(pIndMap,:);
    if size(A,2)>= R.SimAn.minRank
        B = eig(cov(A'));
        C = B/sum(B);
        eRank = sum(cumsum(C)>0.01);
        R.SimAn.minRank = ceil(eRank*2);
        fprintf('effective rank of optbank is %.0f\n',eRank)
    end
    if size(parOptBank,2)> R.SimAn.minRank-1
        if size(parOptBank,2) < (R.SimAn.minRank-1)
            disp('Bank satisfies current eps')
            eps_act = eps_exp;
            cflag = 1; % copula flag (enough samples)
            itry = 0;  % set counter to 0
        else % if the bank is very large than take subset
            disp('Bank is large taking new subset to form eps')
            parOptBank = parBank(:,intersect(1:R.SimAn.minRank,1:size(parBank,2)));
            ACClocbank = computeObjective(R,parOptBank(end,:));
            eps_act = prctile(ACClocbank,50);
            cflag = 1; % copula flag (enough samples)
            itry = 0;  % set counter to 0
        end
    elseif (itry < 1) || (size(parBank,2) < (R.SimAn.minRank-1))
        fprintf('Trying for the %.0f\n time with the current eps \n',itry)
        disp('Trying once more with current eps')
        if isfield(Mfit,'Rho')
            cflag = 1;
        end
        itry = itry + 1;
    elseif itry >= 1
        disp('Recomputing eps from parbank')
        parOptBank = parBank(:,intersect(1:R.SimAn.minRank,1:size(parBank,2)));
        ACClocbank = computeObjective(R,parOptBank(end,:));
        eps_act = prctile(ACClocbank,50);
        cflag = 1;
        itry = 0;
        
        % Stop getting stuck
        %         if (eps_act-eps_prior) == 0
        %             eps_act = eps_act./2;
        %         end
    end
    
    if itry==0
        % Compute expected gradient for next run
        delta_exp = eps_exp-eps_prior;
        fprintf('Expected gradient was %0.2f \n',delta_exp)
        delta_act = eps_act-eps_prior;
        fprintf('Actual gradient was %0.2f \n',delta_act)
        eps_exp = eps_act + delta_act;
        fprintf('Exp-Act gradient was %0.2f \n',delta_exp-delta_act)
        % Save eps history and make actual eps new prior eps
        eps_prior = eps_act;
    end
    eps_rec(ii) = eps_act;
    
    %% Compute Proposal Distribution
    if cflag == 1 && itry == 0 % estimate new copula
        [Mfit,cflag] = postEstCopula(parOptBank,Mfit,pIndMap,pOrg);
        [KL,DKL,R] = KLDiv(R,Mfit,pOrg,m,1);
        Mfit.DKL = DKL;
    elseif cflag == 0 && (size(parBank,2) >= (R.SimAn.minRank-1))% estimate mv Normal Distribution
        % Set Weights
        if size(parOptBank,2)>R.SimAn.minRank
            s = parOptBank(end,:);
            xs = parOptBank(pIndMap,:);
        else
            s = parBank(end,intersect(1:R.SimAn.minRank,1:size(parBank,2)));
            xs = parBank(pIndMap,intersect(1:R.SimAn.minRank,1:size(parBank,2)));
        end
        W = ((s(end,:)-1).^-1);
        W = W;%.^2;
        W = W./sum(W);
        Ws = repmat(W,size(xs,1),1); % added 03/2020 as below wasnt right dim (!)
        Mfit.Mu = wmean(xs,Ws,2);
        Mfit.Sigma = weightedcov(xs',W);
        R.Mfit = Mfit;
        
        %             Mfit.Mu = mean(parBank(pMuMap,intersect(1:1.5*R.SimAn.minRank,1:size(parBank,2))),2);
        %             Mfit.Sigma = cov(parBank(pMuMap,intersect(1:1.5*R.SimAn.minRank,1:size(parBank,2)))');
        
        [KL,DKL,R] = KLDiv(R,Mfit,pOrg,m,0);
        Mfit.DKL = DKL;
    end
    
    %% Draw from Proposal Distribution
    if cflag == 1
        [par,MAP] = postDrawCopula(R,Mfit,pOrg,pIndMap,pSigMap,rep);
        Mfit.MAP = MAP;
        R.Mfit = Mfit;
    elseif cflag == 0
        par = postDrawMVN(R,Mfit,pOrg,pIndMap,pSigMap,rep);
    end
    try
        kldHist(ii) = R.SimAn.scoreweight(2)*R.Mfit.DKL;
        r2Hist(ii) = mean(R.SimAn.scoreweight(2).*r2rep);
    catch
        kldHist(ii) = NaN;
        r2Hist(ii) = NaN;
    end
    saveMkPath([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\klHist_' R.out.tag '_' R.out.dag '.mat'],kldHist)
    parPrec(:,ii+1) = diag(Mfit.Sigma);
    deltaPrec(ii) = mean(diff(parPrec(:,[ii ii+1]),[],2));
    parHist(ii) = averageCell(par);
    saveMkPath([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag  '\' R.out.dag '\parHist_' R.out.tag '_' R.out.dag '.mat'],parHist)
    banksave{ii} = parBank(end,parBank(end,:)>eps_act);
    saveMkPath([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag  '\' R.out.dag '\bankSave_' R.out.tag '_' R.out.dag '.mat'],banksave)
    %%%%%%%%%%%%%%% SAVE PROGRESS, PLOTTING ETC. %%%%%%%%%%%%%%%%%%%%%%%%%%
    if R.plot.flag ==1
        if isfield(R.plot,'outFeatFx')
            %% Plot Data Features Outputs
            try
                R.plot.outFeatFx({R.data.feat_emp},bestfeat,R.data.feat_xscale,R,1,[])
                drawnow; %shg
            catch
                disp('Time series plotting failed!')
            end
        end
        %% Plot parameters changes and tracking of fit
        if isfield(Mfit,'Rho')
            pmean = Mfit.Pfit;
        else
            pmean = p;
        end
        
        figure(2);  clf
        optProgPlot(1:ii,bestr2,pmean,banksave,eps_rec,bestr2,pInd,pSig,R,kldHist,r2Hist)
        drawnow;%shg
        %% Plot example time series
        try
            figure(22)
            plotTimeSeriesGen(xsims_gl_best,1./R.IntP.dt,R.chsim_name,R.condnames)
        catch
             disp('Feature plotting failed!')
        end
    end
    disp({['Current R2: ' num2str(bestr2)];[' Temperature ' num2str(ii) ' K']; R.out.tag; ['Eps ' num2str(eps)]})
    
    %% Save data
    if rem(ii,10) == 0 || ii == 1
        saveSimABCOutputs(R,Mfit,m,parBank)
        if R.plot.save == 1
            saveSimAnFigures(R,ii)
        end
    end
    
    try
        RFLAG = (numel(unique(eps_rec(end-R.SimAn.convIt.eqN:end))) == 1);
    catch
        RFLAG = 0;
    end
    %
    if ~rem(ii,5)
    a= 1;    
    end
    
    if (abs(delta_act) < R.SimAn.convIt.dEps && abs(delta_act)~=0) || RFLAG
        disp('Itry Exceeded: Convergence')
        saveSimABCOutputs(R,Mfit,m,parBank)
        if R.plot.flag == 1
            H(1) = figure(1);
            H(2) = figure(2);
            H(3) = figure(3);
            saveFigure(H,[R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\convergenceFigures'])
        end
        return
    end
    
    ii = ii + 1;
    %     uv = whos;
    %     saveMkPath([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\memDebug_' R.out.tag '_' R.out.dag '.mat'],uv)
    %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%    %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
end
