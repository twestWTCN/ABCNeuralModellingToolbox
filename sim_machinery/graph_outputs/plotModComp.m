function plotModComp_091118(R,permMod)
% addpath('C:\Users\Tim\Documents\MATLAB_ADDONS\violin')
R.plot.confint = 'none';
% cmap = brewermap(7,'Set2');
cmap = linspecer(size(permMod,2));
% cmap(7,:) = cmap(8,:);
% cmap(8,:) = cmap(9,:);
% 
mni = 0;
for mnum = 1:size(permMod,2)
    A = permMod{mnum};
    
    if ~isempty(A)
        r2rep = [A.r2rep{:}];
        r2rep(isnan(r2rep) | isinf(r2rep)) = [];
        
        r2repSave{mnum} = r2rep;
        pd = fitdist(r2rep','Normal');
        x_values = -1:0.01:1;
        y = pdf(pd,x_values);
        figure(1)
        plot(x_values,y,'LineWidth',2)
        
        %         figure(1)
        %         histogram(r2rep,-1:0.1:1);
        hold on
        
        KL(mnum) = sum(A.KL);
        DKL(mnum) = sum(A.DKL);
        pmod(mnum) = sum(r2rep>R.analysis.modEvi.eps)/ size(r2rep,2);
        
        h = figure(10);
        R.plot.cmap = cmap(mnum,:);
        [hl(mnum), hp, dl, flag] = PlotFeatureConfInt_gen060818(R,permMod{mnum},h);
        if ~flag
            mni = mni +1;
            longlab{mni} = sprintf('Model %.f',mnum);
        end
    else
        r2repSave{mnum} = nan(size(r2rep));
    end
    
    shortlab{mnum} = sprintf('M%.f',mnum);
end
hl(end+1) = dl;
longlab{end+1} = 'Data';
figure(30)
violin(r2repSave,'facecolor',cmap,'medc','k:','xlabel',shortlab)
xlabel('Model'); ylabel('NMRSE'); grid on; ylim([-1.5 0.25])

figure(10)
h = findobj(gca,'Type','line');
legend(h([size(h,1):-2:1 1]),longlab)

figure(2)
subplot(2,1,1)
for i = 1:size(permMod,2)
    b = bar(i,pmod(i)); hold on
    b.FaceColor = cmap(i,:);
end
a = gca; a.XTick = 1:size(permMod,2); grid on
a.XTickLabel = shortlab; 
xlabel('Model'); ylabel('P(M|D)')
subplot(2,1,2)
for i = 1:size(permMod,2)
    b = bar(i,KL(i)); hold on
    b.FaceColor = cmap(i,:);
end
a = gca; a.XTick = 1:size(permMod,2);
a.XTickLabel = shortlab; 
grid on
xlabel('Model'); ylabel('KL Divergence')

% subplot(3,1,3)
% bar(DKL)
% xlabel('Model'); ylabel('Joint KL Divergence')