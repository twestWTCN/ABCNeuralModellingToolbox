function plotModComp_310520(R,cmap)
% addpath('C:\Users\Tim\Documents\MATLAB_ADDONS\violin')
% R.plot.confint = 'none';
if nargin<2
    cmap = brewermap(R.modcomp.modN,'Spectral');
    % cmap = linspecer(R.modcomp.modN);
end
%% First get probabilities so you can compute model space epsilon
for modID = 1:numel(R.modcomp.modN)
    R.out.dag = sprintf([R.out.tag '_M%.0f'],modID);
    load([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\modeProbs_' R.out.tag '_'  R.out.dag '.mat'])
    permMod = varo; %i.e. permMod
    if ~isempty(permMod)
        ACCrep{modID} = permMod.ACCrep;
    end
end

% This gets the joint space epsilon that you can use to compute exceedence
% probability
prct = 50;
ACCbankcat = horzcat(ACCrep{:});
R.modcomp.modEvi.epspop = prctile(ACCbankcat,prct); % threshold becomes median of model fits

p = 0; % plot counter
mni = 0; % Model counter
for modID = 1:numel(R.modcomp.modN)
    shortlab{modID} = sprintf('M%.f',R.modcomp.modN(modID)); % Make model label

    % Load in the precompute model iterations
    R.out.dag = sprintf([R.out.tag '_M%.0f'],R.modcomp.modN(modID));
    load([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\modeProbs_' R.out.tag '_'  R.out.dag '.mat'])
    permMod = varo; %i.e. permMod

    load([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\R_' R.out.tag '_' R.out.dag  '.mat'])
    tmp = varo;
    R.data = tmp.data;
    if ~isempty(permMod)
        % save the model accuracies
        ACCrep = permMod.ACCrep;

        % now get the exceedence probability
        pmod(modID) =sum(ACCrep>R.modcomp.modEvi.epspop) / size(ACCrep,2);
        % parameter divergences
        KL(modID) = sum(permMod.KL(~isnan(permMod.KL))); % sum across all (marginal distributions)
        DKL(modID) = permMod.DKL; % total joint space KL divergence
        MAP(modID) = permMod.MAP;
        ACC{modID} = permMod.ACCrep;

        %% Plot Data Features with Bayesian confidence intervals
        h(1,1) = figure(10);
        h(2,1) = figure(20);
        flag = 0;

        if ismember(modID,R.modcompplot.NPDsel)
            p = p +1;
            [hl(p), hp, dl, flag] = PlotFeatureConfInt_gen170620(R,permMod,h, cmap(modID,:));
        end
        % hl(modID) = plot(1,1);
        if ~flag
            mni = mni +1;
        end
    else
        r2repSave{modID} = nan(size(ACCrep));
    end
end
% Save the model parameter average
save([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.tag '_model_parameter_averages'],'parMean')
set(gcf,'Position',[680   112   976   893])


%% Now Plot Results of Model Comparison
figure(2)
subplot(4,1,1)
violin(r2repSave,'facecolor',cmap,'medc','k:','xlabel',shortlab); %,...
%'bw',[0.025 0.1 0.1]); % 0.025 0.2 0.2 0.025 0.025 0.025 0.025 0.2 0.2]);
hold on
plot([0 numel(R.modcomp.modN)+1],[R.modcomp.modEvi.epspop R.modcomp.modEvi.epspop],'k--')
xlabel('Model'); ylabel('NMRSE'); grid on;
% ylim([-1 0.4])
a = gca; a.XTick = 1:numel(R.modcomp.modN);
a.XTickLabel = shortlab;

figure(2)
subplot(4,1,1)
h = findobj(gca,'Type','line');
% legend(hl,{longlab{[R.modcompplot.NPDsel end]}})

subplot(4,1,2)
% TlnK = 2.*log(max(pmod)./pmod);
%  TlnK = log10(pmod); % smallest (closest to zero) is the best
TlnK = -log10(1-pmod); % largest is the best

for i = 1:numel(R.modcomp.modN)
    %     b = bar(i,-log10(1-pmod(i))); hold on
    b = bar(i,TlnK(i)); hold on
    b.FaceColor = cmap(i,:);
end
a = gca; a.XTick = 1:numel(R.modcomp.modN); grid on
a.XTickLabel = shortlab;
xlabel('Model'); ylabel('-log_{10} P(M|D)')
xlim([0.5 numel(R.modcomp.modN)+0.5])

subplot(4,1,3)
for i = 1:numel(R.modcomp.modN)
    b = bar(i,DKL(i)); hold on
    b.FaceColor = cmap(i,:);
end
a = gca; a.XTick = 1:numel(R.modcomp.modN);
a.XTickLabel = shortlab;
grid on
xlabel('Model'); ylabel('KL Divergence')
set(gcf,'Position',[277   109   385   895])
xlim([0.5 numel(R.modcomp.modN)+0.5])
% subplot(3,1,3)
% bar(DKL)
% xlabel('Model'); ylabel('Joint KL Divergence')

subplot(4,1,4)
for i = 1:numel(R.modcomp.modN)
    b = bar(i,ACC(i)); hold on
    b.FaceColor = cmap(i,:);
end
a = gca; a.XTick = 1:numel(R.modcomp.modN);
a.XTickLabel = shortlab;
grid on
xlabel('Model'); ylabel('KL Divergence')
set(gcf,'Position',[277   109   385   895])
xlim([0.5 numel(R.modcomp.modN)+0.5])



%% SCRIPT GRAVE
% Adjust the acceptance threshold if any models have no rejections
% exc = ones(1,numel(R.modcomp.modN));
% while any(exc==1)
%     r2bankcat = horzcat(r2bank{:});
%     R.modcomp.modEvi.epspop = prctile(r2bankcat,prct); % threshold becomes median of model fits
%     for modID = 1:numel(R.modcomp.modN)
%         exc(modID) = sum(r2bank{modID}>R.modcomp.modEvi.epspop)/size(r2bank{modID},2);
%     end
%     prct = prct+1;
% end
