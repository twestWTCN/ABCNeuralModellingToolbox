function [R, parBank, Mfit] = SimAn_ABC_SMC_140726(R, p, m, parBank)
%%%% APPROXIMATE BAYESIAN COMPUTATION — SEQUENTIAL MONTE CARLO
%%%% HIGH DIMENSIONAL DYNAMICAL MODELS
% ---- 14/07/26 --------------------------------------------------
% Proper ABC-SMC built on top of SimAn_ABC_140726:
%
%  (1) Fixed population of Npop accepted particles per generation.
%      parBank is kept as a full history for output/plotting only;
%      only the current generation's particles are used to fit the
%      next proposal.
%
%  (2) Importance weights  w_i = pi(theta_i) / q_t(theta_i)
%      where q_t is the fitted t-copula with KDE marginals.
%      All weights are tracked in log-space throughout.
%      For the MVN fallback path, the MVN density is used as q_t.
%
%  (3) Effective Sample Size is tracked each generation.
%      ESS = 1 / sum(w_i^2).  Low ESS signals weight degeneracy.
%
%  (4) GA sample selection uses IS-weighted statistics instead of
%      raw score-based weights.
%
%  (5) Epsilon is computed as the IS-weighted 25th percentile of
%      generation scores, making it scale-independent.
%
%  (6) All helper functions are local (self-contained file).
%
% Timothy West — UCL (2018); SMC revision 2026
%%%% ----------------------------------------------------------------
%%%  %%%  %%%  %%%  %%%  %%%  %%%  %%%  %%%  %%%  %%%  %%%  %%%  %%%
GPool = gcp;
warning('off', 'MATLAB:MKDIR:DirectoryExists');
R = ABCForwardCompatibility(R);
ABCGraphicsDefaults(R)

%% Defaults
if nargin < 4,                        parBank = [];             end
if ~isfield(R.plot,  'flag'),         R.plot.flag       = 1;    end
if ~isfield(R.plot,  'save'),         R.plot.save       = 0;    end
if isempty(m),                        m.m               = 1;    end
if ~isfield(R.plot,  'updateflag'),   R.plot.updateflag = 0;    end
if ~isfield(R.SimAn, 'minRankLambda'),R.SimAn.minRankLambda = 3;end
if ~isfield(R.SimAn, 'RealzRep'),     R.SimAn.RealzRep  = 1;   end
if ~isfield(R.SimAn, 'sigmaAlpha'),   R.SimAn.sigmaAlpha= 0.01; end
if ~isfield(R.SimAn, 'rhoAlpha'),     R.SimAn.rhoAlpha  = 0.1;  end

pOrg = p;

eps_prior = -200;
eps_exp   = -12;
eps_act   = eps_prior;
delta_act = 0.05;
delta_actPerc = 0;

[pInd,~,pSig,pIndMap,pMuMap,pSigMap] = parOptInds_110817(R, p, m.m);
R.SimAn.minRank = ceil(size(pIndMap,1) * R.SimAn.minRankLambda);
minRankFloor    = R.SimAn.minRank;

if ~isfield(R.SimAn,'Npop') || isempty(R.SimAn.Npop)
    R.SimAn.Npop = minRankFloor;   % fixed generation population size
end

%% Initialise proposal distribution
if isfield(R, 'Mfit')
    Mfit           = R.Mfit;
    Mfit.prior     = Mfit;
    Mfit.DKL       = 0;
    cflag_proposal = isfield(Mfit, 'Rho');
else
    ptmp           = spm_vec(p);
    Mfit.Mu        = ptmp(pMuMap);
    Mfit.Sigma     = diag(ptmp(pSigMap) .* R.SimAn.jitter);
    Mfit.prior     = Mfit;
    Mfit.DKL       = 0;
    cflag_proposal = 0;
end

% proposalMfit is the Mfit that generated the CURRENT proposals par.
% It is saved before updating Mfit so IS weights can be computed
% against the correct q_t.
proposalMfit = Mfit;

R.Mfit = Mfit;
parPrec(:,1) = diag(Mfit.Sigma);
cflag = 0;
ii    = 1;

% Draw initial proposals
rep    = R.SimAn.rep(1);
repBig = max(rep, R.SimAn.Npop * 4);   % large enough that circular
                                        % reuse is rare
if cflag_proposal
    [par, ~] = postDrawCopulaPerm(R, proposalMfit, pOrg, pIndMap, pSigMap, repBig, 0);
else
    par = postDrawMVN(R, proposalMfit, pOrg, pIndMap, pSigMap, repBig);
end

%%%  %%%  %%%  %%%  %%%  %%%  %%%  %%%  %%%  %%%  %%%  %%%  %%%  %%%
%% Main Annealing Loop
while ii <= R.SimAn.searchMax

    %% Accumulate generation population
    genBank  = [];      % (nParams+1) x Npop  — current generation
    genLogW  = [];      % 1 x Npop            — log IS weights
    featbank = {}; ACCbank = []; samppar = {};
    ji       = 0;
    parnum   = 4 * GPool.NumWorkers;

    while size(genBank, 2) < R.SimAn.Npop

        % Circular index into proposal batch
        parl_base = mod(ji * parnum, repBig);

        clear r2rep ACCrep par_rep feat_sim_rep
        parfor jj = 1:parnum
            parl = mod(parl_base + jj - 1, repBig) + 1;
            pnew = par{parl};
            r2 = []; feat_sim = [];
            for RzRep = 1:R.SimAn.RealzRep
                [r2(RzRep), pnew, feat_sim] = computeSimData_160620(R, m, [], pnew, 0, 0);
            end
            r2 = mean(r2);
            [ACC, R2w]       = computeObjective(R, r2);
            r2rep(jj)        = R2w;
            ACCrep(jj)       = ACC;
            par_rep{jj}      = pnew;
            feat_sim_rep{jj} = feat_sim;
        end

        if rem(ji, 4) == 1
            fprintf('Batch %d | gen %d | pop %d/%d\n', ji, ii, size(genBank,2), R.SimAn.Npop)
        end

        % Collect valid results into a matrix
        r2loop = r2rep;
        r2loop(r2loop == 1)     = -inf;
        r2loop(isnan(r2loop))   = -inf;
        r2loop(~isreal(r2loop)) = -inf;
        r2loop(isinf(r2loop))   = -inf;

        newParMat = [];
        for i = 1:parnum
            if ~isinf(r2loop(i))
                newParMat = [newParMat, [full(spm_vec(par_rep{i})); r2loop(i)]]; %#ok<AGROW>
            end
        end

        % --- Batch IS weight computation (outside parfor) ---
        if ~isempty(newParMat)
            if cflag_proposal   % t-copula proposal
                log_q = computeCopulaDensity(newParMat, proposalMfit, pIndMap);
            else                % MVN proposal
                log_q = computeMVNDensity(newParMat, proposalMfit, pIndMap);
            end
            log_pi = computePriorDensity(newParMat, pIndMap, proposalMfit.prior);
            log_w  = log_pi - log_q;

            % Cap extreme weights at 95th-pct + log(20) to prevent
            % single particles dominating
            if numel(log_w) > 1
                log_w = min(log_w, quantile(log_w, 0.95) + log(20));
            end

            genBank = [genBank,  newParMat]; %#ok<AGROW>
            genLogW = [genLogW,  log_w];     %#ok<AGROW>
        end

        % Append to history bank and trim
        parBank = [parBank, newParMat]; %#ok<AGROW>
        if ~isempty(parBank)
            [~, V] = sort(parBank(end,:), 'descend');
            if size(parBank,2) > 2^10
                parBank = parBank(:, V(1:2^10));
            else
                parBank = parBank(:, V);
            end
        end

        featbank{ji+1}  = feat_sim_rep; %#ok<AGROW>
        ACCbank(:,ji+1) = ACCrep;       %#ok<AGROW>
        samppar{ji+1}   = par_rep;      %#ok<AGROW>
        ji = ji + 1;

        % Refresh proposal batch when circular window wraps
        if mod(ji * parnum, repBig) < parnum
            if cflag_proposal
                [par, ~] = postDrawCopulaPerm(R, proposalMfit, pOrg, pIndMap, pSigMap, repBig, 0);
            else
                par = postDrawMVN(R, proposalMfit, pOrg, pIndMap, pSigMap, repBig);
            end
        end

    end  % inner while

    %% Trim generation to Npop (keep highest-scoring particles)
    if size(genBank, 2) > R.SimAn.Npop
        [~, V]  = sort(genBank(end,:), 'descend');
        V       = V(1:R.SimAn.Npop);
        genBank = genBank(:, V);
        genLogW = genLogW(V);
    end

    %% Normalise IS weights and compute ESS
    log_Z  = logSumExp(genLogW);
    genW   = exp(genLogW - log_Z);        % normalised, sum to 1
    ESS_ii = 1 / sum(genW .^ 2);
    fprintf('Gen %d: ESS = %.1f / %d (%.0f%%)\n', ii, ESS_ii, R.SimAn.Npop, 100*ESS_ii/R.SimAn.Npop)

    %% Find best draws for plotting
    bestfeat = [];
    [~, Ib] = maxk(ACCbank(:), 12);
    [jj_best, ji_best] = ind2sub(size(ACCbank), Ib);
    pnew = samppar{ji_best(1)}{jj_best(1)};
    [~,~,~,~,xsims_gl_best] = computeSimData_160620(R, m, [], pnew, 0, 0);
    for L = 1:numel(Ib)
        for j = 1:numel(featbank{ji_best(L)}{jj_best(L)})
            bestfeat{L}{j} = featbank{ji_best(L)}{jj_best(L)}{j}; %#ok<AGROW>
        end
    end
    bestr2(ii) = genBank(end, 1); %#ok<AGROW>

    %% IS-weighted 25th-percentile epsilon
    [sorted_s, sort_idx] = sort(genBank(end,:), 'ascend');
    cum_w   = cumsum(genW(sort_idx));
    hit     = find(cum_w >= 0.25, 1, 'first');
    if isempty(hit), hit = 1; end
    eps_act = sorted_s(hit);

    delta_exp     = eps_exp - eps_prior;
    delta_act     = eps_act - eps_prior;
    delta_actPerc = 100 * ((eps_act - eps_prior) / (abs(eps_prior) + eps));
    eps_exp       = eps_act + delta_act;
    fprintf('Eps: prior=%.3f  act=%.3f  exp=%.3f  pct-change=%.1f%%\n', ...
        eps_prior, eps_act, eps_exp, delta_actPerc)
    eps_prior = eps_act;
    eps_rec(ii) = eps_act; %#ok<AGROW>

    %% Effective rank and minRank floor
    A = genBank(pIndMap,:);
    B = eig(cov(A'));
    C = B / sum(B);
    eRank = sum(cumsum(C) > 0.01);
    R.SimAn.minRank = max(ceil(eRank * R.SimAn.minRankLambda), minRankFloor);
    fprintf('Effective rank: %.0f  minRank: %.0f (floor %.0f)\n', eRank, R.SimAn.minRank, minRankFloor)

    %% GA selection (IS-weighted)
    Mfit.W = genW;
    optIdx = find(selectIndicesGASMC(R, genBank, pIndMap, pOrg, genW));
    if numel(optIdx) < R.SimAn.minRank
        % GA failed to find valid subset — fall back to top-weight particles
        [~, sortW] = sort(genW, 'descend');
        optIdx = sortW(1:min(R.SimAn.minRank, numel(sortW)));
    end
    parOptBank = genBank(:, optIdx);
    optW       = genW(optIdx);
    optW       = optW / sum(optW);           % renormalise after selection
    Mfit.W     = optW;

    %% Fit proposal distribution
    Mfit.rhoAlpha = R.SimAn.rhoAlpha;
    [Mfit, cflag] = postEstCopulaSMC(parOptBank, Mfit, pIndMap, pOrg);

    Mfit.Sigma = floorSigma(Mfit.Sigma, Mfit.prior.Sigma, R.SimAn.sigmaAlpha);
    Mfit.Sigma = shrinkCorr(Mfit.Sigma, R.SimAn.rhoAlpha);

    if cflag   % copula path
        [KL, DKL, R] = KLDiv(R, Mfit, pOrg, m, 1); %#ok<ASGLU>
    else       % MVN fallback (copulafit failed)
        R.Mfit = Mfit;
        [~, DKL, R] = KLDiv(R, Mfit, pOrg, m, 0);
    end
    Mfit.DKL = DKL;
    R.Mfit   = Mfit;

    %% Draw proposals for the NEXT generation
    % Save proposalMfit so IS weights can be computed against q_t next gen
    proposalMfit   = Mfit;
    cflag_proposal = cflag;

    if cflag
        [par, MAP] = postDrawCopulaPerm(R, proposalMfit, pOrg, pIndMap, pSigMap, repBig, 0);
        Mfit.MAP   = MAP;
    else
        par = postDrawMVN(R, proposalMfit, pOrg, pIndMap, pSigMap, repBig);
    end

    %% Diagnostics
    try
        kldHist(ii) = R.SimAn.scoreweight(2) * Mfit.DKL; %#ok<AGROW>
        r2Hist(ii)  = genW * genBank(end,:)';             %#ok<AGROW>  IS-weighted mean score
        ESSHist(ii) = ESS_ii;                             %#ok<AGROW>
    catch
        kldHist(ii) = NaN; r2Hist(ii) = NaN; ESSHist(ii) = NaN; %#ok<AGROW>
    end

    saveMkPath([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\klHist_'  R.out.tag '_' R.out.dag '.mat'], kldHist)
    saveMkPath([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\ESSHist_' R.out.tag '_' R.out.dag '.mat'], ESSHist)

    parPrec(:,ii+1) = diag(Mfit.Sigma); %#ok<AGROW>
    parHist(ii)     = averageCell(par);  %#ok<AGROW>
    saveMkPath([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\parHist_' R.out.tag '_' R.out.dag '.mat'], parHist)

    banksave{ii} = parBank(end, parBank(end,:) > eps_act); %#ok<AGROW>
    saveMkPath([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\bankSave_' R.out.tag '_' R.out.dag '.mat'], banksave)

    %% Plotting
    if R.plot.flag == 1
        if isfield(R.plot, 'outFeatFx')
            try
                set(groot, 'CurrentFigure', 11); clf
                R.plot.outFeatFx({R.data.feat_emp}, bestfeat, R.data.feat_xscale, R, 1, [])
                drawnow
            catch; disp('Feature plotting failed'); end
        end
        if isfield(Mfit, 'Rho'), pmean = Mfit.Pfit; else, pmean = p; end
        try
            set(groot, 'CurrentFigure', 100); clf
            optProgPlot(1:ii, bestr2, pmean, banksave, eps_rec, bestr2, pInd, pSig, R, kldHist, r2Hist)
            drawnow
        end
        try
            set(groot, 'CurrentFigure', 22)
            plotTimeSeriesGen(xsims_gl_best, 1./R.IntP.dt, R.chsim_name, R.condnames)
        catch; disp('Time series plotting failed'); end
    end

    disp({['IS-wtd R2: ' num2str(r2Hist(ii))]; ...
          ['ESS: '       num2str(ESS_ii, '%.1f') '/' num2str(R.SimAn.Npop)]; ...
          [' Gen '       num2str(ii)]; R.out.tag; R.out.dag; ...
          ['Eps change ' num2str(delta_actPerc, '%.1f') '%']})

    if rem(ii,1) == 0 || ii == 1
        saveSimABCOutputs(R, Mfit, m, parBank)
        if R.plot.save == 1, saveSimAnFigures(R, ii); end
    end

    %% Convergence check (fractional epsilon gradient — unitless)
    try
        RFLAG    = (numel(unique(eps_rec(end-R.SimAn.convIt.eqN:end))) == 1);
        epsWin   = eps_rec(end-R.SimAn.convIt.eqN:end);
        meanGrad = mean(diff(epsWin)) / (abs(mean(epsWin)) + eps);
    catch
        meanGrad = eps_act;
        RFLAG    = 0;
    end

    if isfield(R.plot, 'updateperiod')
        if ~rem(ii, R.plot.updateperiod)
            R.plot.flag = 1; ABCGraphicsDefaults
        else
            R.plot.flag = 0;
        end
    end

    if (abs(meanGrad) < R.SimAn.convIt.dEps && abs(delta_act) ~= 0) || RFLAG || ii > R.SimAn.convIt.MaxIt
        disp('Convergence criterion met')
        saveSimABCOutputs(R, Mfit, m, parBank)
        if R.plot.save == 1
            H = get(groot, 'Children');
            saveFigure(H, [R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\' R.out.dag '\convergenceFigures'])
        end
        return
    end

    ii = ii + 1;
    %%%  %%%  %%%  %%%  END OF GENERATION  %%%  %%%  %%%  %%%  %%%  %%%
end


%%%  %%%  %%%  %%%  %%%  LOCAL FUNCTIONS  %%%  %%%  %%%  %%%  %%%  %%%

%% -----------------------------------------------------------------------
function [Mfit, cflag] = postEstCopulaSMC(parOptBank, Mfit, pIndMap, pOrg)
% POSTESTCOPULASMC  Fit t-copula using IS weights (Mfit.W).
% Falls back to score-based weights if Mfit.W is absent (backward compat).
disp('Forming copula (IS-weighted)...')
clear copU xf

% Select weight source
if isfield(Mfit,'W') && numel(Mfit.W) == size(parOptBank,2)
    W = Mfit.W(:)';
    W = W / sum(W);
else
    W = (parOptBank(end,:) - 1).^-1;
    W = min(W, quantile(W, 0.95));
    W = W ./ sum(W);
end

% Per-parameter bandwidth floor from prior covariance
priorStd  = sqrt(abs(diag(Mfit.prior.Sigma)));
bwidFloor = 0.1 * priorStd';

for i = 1:size(pIndMap, 1)
    x       = parOptBank(pIndMap(i), :);
    bwid(i) = KSDensityCVWidth(x, x, W, [-1 1], 25, 'cdf');
    bwid(i) = max(bwid(i), bwidFloor(i));
    copU(i,:) = ksdensity(x, x, 'function','cdf', 'Weights', W, 'width', bwid(i));
    xf(i,:) = x;
end

cflag = 0;
try
    [Rho, nu] = copulafit('t', copU', 'Method','ApproximateML');
    ra   = Mfit.rhoAlpha;
    Rho  = (1-ra)*Rho + ra*eye(size(Rho));
    Rho  = (Rho + Rho') / 2;

    Mfit.xf   = xf;
    Mfit.ks   = copU;
    Mfit.nu   = nu;
    Mfit.bwid = bwid;
    Mfit.tbr2 = parOptBank(end, 1);
    % IS-weighted mean as point estimate
    Mfit.Pfit  = spm_unvec(parOptBank(1:end-1,:) * W(:), pOrg);
    Mfit.BPfit = spm_unvec(parOptBank(1:end-1, 1), pOrg);
    Mfit.Rho   = Rho;
    cflag      = 1;
catch
    disp('Copula estimation failed (rank deficient?) — falling back to MVN')
end

% IS-weighted Mu and Sigma (used for KLD and pSigMap slots)
xs         = parOptBank(pIndMap,:);
Ws         = repmat(W, size(xs,1), 1);
Mfit.Mu    = wmean(xs, Ws, 2);
Mfit.Sigma = weightedcov(xs', W);


%% -----------------------------------------------------------------------
function optimalIndices = selectIndicesGASMC(R, parBank, pIndMap, pOrg, W)
% SELECTINDICESGASMC  GA sample selection using IS weights.
options = optimoptions('ga', ...
    'PopulationSize',   100, ...
    'MaxGenerations',   150, ...
    'EliteCount',       2,   ...
    'CrossoverFraction',0.8, ...
    'MutationFcn',      @mutationuniform, ...
    'Display',          'none');
numSamples = size(parBank, 2);
IntCon     = 1:numSamples;
disp('Using GA to select proposal samples (IS-weighted)')
fitFcn = @(idx) sampleFitnessSMC(idx, R, parBank, pIndMap, pOrg, W);
[optimalIndices, ~] = ga(fitFcn, numSamples, [],[],[],[], ...
    zeros(1,numSamples), ones(1,numSamples), [], IntCon, options);


%% -----------------------------------------------------------------------
function fitness = sampleFitnessSMC(indices, R, parBank, pIndMap, pOrg, W)
% SAMPLEFITNESSSMC  Fitness using IS-weighted mean score and KLD.
selected = find(indices);
if numel(selected) < R.SimAn.minRank
    fitness = inf;
    return
end
parOptBank = parBank(:, selected);
Wsel       = W(selected);
Wsel       = Wsel / (sum(Wsel) + eps);

% IS-weighted mean score
meanError = Wsel * parOptBank(end,:)';

% IS-weighted Mu and Sigma for KLD
xs              = parOptBank(pIndMap,:);
Ws              = repmat(Wsel, size(xs,1), 1);
Mfit_tmp.prior  = R.Mfit.prior;
Mfit_tmp.Mu     = wmean(xs, Ws, 2);
Mfit_tmp.Sigma  = weightedcov(xs', Wsel);
[~, DKL, ~]    = KLDiv(R, Mfit_tmp, pOrg, [], 0);

fitness = -(R.SimAn.scoreweight(1)*meanError - R.SimAn.scoreweight(2)*DKL);


%% -----------------------------------------------------------------------
function log_q = computeCopulaDensity(parMat, Mfit, pIndMap)
% COMPUTECOPULADENSITY  Evaluate log q(theta) for the t-copula proposal.
%
%   log q(theta) = log c(F_1(t_1),...,F_d(t_d); Rho, nu)
%                + sum_k log f_k(t_k)
%
%   where c is the t-copula density, F_k the marginal empirical CDF
%   (KDE-based), and f_k the marginal KDE PDF.

xf        = Mfit.xf;
bwid      = Mfit.bwid;
Rho       = Mfit.Rho;
nu        = Mfit.nu;
theta_opt = parMat(pIndMap, :);   % d x nSamples
nS        = size(theta_opt, 2);
d         = size(xf, 1);

U             = zeros(nS, d);
log_marginals = zeros(nS, d);

for k = 1:d
    Uk = ksdensity(xf(k,:), theta_opt(k,:)', 'function','cdf', 'width', bwid(k));
    fk = ksdensity(xf(k,:), theta_opt(k,:)', 'function','pdf', 'width', bwid(k));
    U(:,k)             = Uk;
    log_marginals(:,k) = log(max(fk, 1e-300));
end

% Clip U strictly inside (0,1) to keep copulapdf finite
U = min(max(U, 1e-6), 1-1e-6);

c     = copulapdf('t', U, Rho, nu);
log_q = (log(max(c, 1e-300)) + sum(log_marginals, 2))';   % 1 x nS


%% -----------------------------------------------------------------------
function log_q = computeMVNDensity(parMat, Mfit, pIndMap)
% COMPUTEMVNDENSITY  Log MVN proposal density (used before copula fitted).
theta_opt = parMat(pIndMap, :)';                            % nS x d
log_q     = log(mvnpdf(theta_opt, Mfit.Mu(:)', Mfit.Sigma) + 1e-300)';


%% -----------------------------------------------------------------------
function log_pi = computePriorDensity(parMat, pIndMap, prior)
% COMPUTEPRIORDENSITY  Log Gaussian prior density.
theta_opt = parMat(pIndMap, :)';                            % nS x d
log_pi    = log(mvnpdf(theta_opt, prior.Mu(:)', prior.Sigma) + 1e-300)';


%% -----------------------------------------------------------------------
function lse = logSumExp(x)
% LOGSUMEXP  Numerically stable log(sum(exp(x))).
m   = max(x);
lse = m + log(sum(exp(x - m)));


%% -----------------------------------------------------------------------
function Sigma = floorSigma(Sigma, priorSigma, alpha)
% FLOORSIGMA  Hard eigenvalue floor at alpha * min(prior eigenvalue).
priorEigs = eig(priorSigma);
eigFloor  = alpha * min(priorEigs(priorEigs > 0));
[V, D]    = eig(Sigma);
d         = max(diag(D), eigFloor);
Sigma     = (V * diag(d) * V' + (V * diag(d) * V')') / 2;


%% -----------------------------------------------------------------------
function Sigma = shrinkCorr(Sigma, alpha)
% SHRINKCORR  Shrink off-diagonal correlations toward zero by factor alpha.
% Marginal variances (diagonal) are unchanged.
stds  = sqrt(diag(Sigma));
Corr  = Sigma ./ (stds * stds');
Corr  = (1-alpha)*Corr + alpha*eye(size(Corr));
Sigma = Corr .* (stds * stds');
Sigma = (Sigma + Sigma') / 2;
