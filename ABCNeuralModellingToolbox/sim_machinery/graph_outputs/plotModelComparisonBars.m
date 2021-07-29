function plotModelComparisonBars(R,compStruc,shortlab,type)
% TlnK = 2.*log(max(pmod)./pmod);
if type == 1 % model probs
    F = -log10(1-compStruc.pmod);
    yl = '-log_{10} P(M|D)';
elseif type == 2 % KL divergences
    F = compStruc.KL;
    yl = 'KL Divergence';
    elseif type == 3 % ACS
    F = compStruc.ACS;
    yl = 'ACS Combined Score';
end

for i = 1:numel(R.modcomp.modN)
    %     b = bar(i,-log10(1-pmod(i))); hold on
    b = bar(i,F(i)); hold on
    b.FaceColor = R.plot.cmap(i,:);
end
a = gca; a.XTick = 1:numel(R.modcomp.modN); grid on
a.XTickLabel = shortlab;
xlabel('Model'); 
ylabel(yl)