function [permMod] = modelBetaAlignedAnaly(x,m,p,R)
% load([R.rootn 'outputs\' R.out.tag '\parBank_' R.out.tag '_' d '.mat'])
parOptBank = R.parOptBank;
% figure
% hist(parOptBank(end,:),[-1:.1:1]); xlim([-1 1])
eps = R.analysis.modEvi.eps;

% parOptBank = parOptBank(:,parOptBank(end,:)>eps);
%% Compute KL Divergence
[KL DKL] = KLDiv(R,p,m,parOptBank)
N = R.analysis.modEvi.N;

%% Resample parameters
% Compute indices of optimised parameter
pInd = parOptInds_110817(R,p,m.m,2); % in structure form
pIndMap = spm_vec(pInd); % in flat form
% pIndMap (71 x 1)
% xf (71 x size(parOptBank,2) )
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
    par{i} = spm_unvec(mean(base,2),p);
end
 plotDistChange_KS(R.Mfit.Rho,R.Mfit.nu,xf,p,pInd,R,1)
% if isempty(gcp)
%     parpool
% end
gcp
ppm = ParforProgMon('Model Probability Calculation',N);
%%Plot Example
figure(5)
pnew = par{2};
%% Simulate New Data
u = innovate_timeseries(R,m);
u{1} = u{1}.*sqrt(R.IntP.dt);
[xsims,tvec,wflag] = R.IntP.intFx(R,m.x,u,pnew,m);
if wflag ==0
    % Run Observer function
    if isfield(R.obs,'obsFx')
        [xsims R] = R.obs.obsFx(xsims,m,pnew,R);
    end
    xsims{1} = xsims{1} + linspace(m.m,0,m.m)';
    plot(R.IntP.tvec_obs(1:end),xsims{1}(:,2:end))
    legend(R.chsim_name); xlabel('Time (s)'); ylabel('Amplitude')
    set(gcf,'Position',[705         678        1210         420]);
    xlim([4 5])
    %%
    close all
end
for jj = 1:N
    %     ppm.increment();
    pnew = par{jj};
    %% Simulate New Data
    u = innovate_timeseries(R,m);
    u{1} = u{1}.*sqrt(R.IntP.dt);
    [xsims,tvec,wflag] = R.IntP.intFx(R,m.x,u,pnew,m);
    if wflag == 0
        % Run Observer function
        if isfield(R.obs,'obsFx')
            xsims = R.obs.obsFx(xsims,m,pnew,R);
        end
        % Run Data Transform
        if isfield(R.obs,'transFx')
            [~,feat_sim] = R.obs.transFx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,R.obs.SimOrd,R);
        else
            feat_sim = xsims; % else take raw time series
        end
        % Compare Pseudodata with Real
        r2mean  = R.IntP.compFx(R,feat_sim);
    else
        r2mean = -inf;
        feat_sim = NaN;
        disp('Simulation error!')
    end
    %     R.plot.outFeatFx({},{feat_sim},R.data.feat_xscale,R,1)
    r2rep{jj} = r2mean;
    par_rep{jj} = pnew;
    feat_rep{jj} = feat_sim;
    %     disp(jj); %
    ppm.increment();
end
permMod.r2rep = r2rep;
permMod.par_rep = par_rep;
permMod.feat_rep = feat_rep;
permMod.DKL = DKL;
permMod.KL = KL;
% mkdir([R.rootn 'outputs\' R.out.tag '2\'])
% save([R.rootn 'outputs\' R.out.tag '2\permMod_' R.out.tag '_' d '.mat'],'permMod')
% load([R.rootn 'outputs\' R.out.tag '2\permMod_' R.out.tag '_' d '.mat'],'permMod')

figure
r2bank = [permMod.r2rep{:}];
[h r] = hist(r2bank,50); %D is your data and 140 is number of bins.
h = h/sum(h); % normalize to unit length. Sum of h now will be 1.
bar(h, 'DisplayName', 'Model NRMSE');
xD = r(2:2:end);
xL = 2:2:length(r); % list of indices
set(gca,'XTick',xL)
set(gca,'XTickLabel',strsplit(num2str(xD,2),' '))

legend('show');
ylabel('P(D-D*)'); xlabel('D-D*');
hold on
Yval = get(gca,'YLim')

tmp = abs(xD-eps);
[idx idx] = min(tmp); %index of closest value
epsm = xL(xD==xD(idx)); %closest value

plot([epsm epsm],Yval,'B--','linewidth',3)

Pmod = numel(r2bank(r2bank>eps))/R.analysis.modEvi.N;
annotation(gcf,'textbox',...
    [0.28 0.81 0.19 0.09],...
    'String',{sprintf('eps = %.2f',eps),sprintf('P(m|D) = %.2f',Pmod)},...
    'HorizontalAlignment','right',...
    'FitBoxToText','off',...
    'LineStyle','none');
set(gcf,'Position',[680 437 1070 541])
% saveallfiguresFIL_n([R.rootn 'outputs\' R.out.tag '\modelEvidence.jpg'],'-jpg',1,'-r200',1);

