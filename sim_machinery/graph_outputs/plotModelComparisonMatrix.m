function h = plotModelComparisonMatrix(R,compStruc,shortlab,type)
% TlnK = 2.*log(max(pmod)./pmod);
if type == 1 % model probs
    F = compStruc.pmod; %-log10(1-compStruc.pmod);
    F = reshape(F,sqrt(numel(F)),sqrt(numel(F)));
    Flab = F;
    F = -log10(1-F);
    yl = 'Posterior Model Probability';
elseif type == 2 % KL divergences
    F = compStruc.DKL;
    F = reshape(F,sqrt(numel(F)),sqrt(numel(F)));
    Flab = F;
    F = -F;
    yl = 'KL Divergence';
elseif type == 3 % Score
    F = compStruc.ACS;
    F = reshape(F,sqrt(numel(F)),sqrt(numel(F)));
%     F = -F;
    Flab = F;
    yl = 'Accuracy-Divergence';
end

h = imagesc(F);
[x,y] = meshgrid(1:size(F,1),1:size(F,2));
text(x(:),y(:),num2str(Flab(:),'%0.3f'),'HorizontalAlignment','center')
cmap = brewermap(128,'*RdYlBu');
colormap(cmap);
a = gca;
% set(a,'YDir','rev') %'XTick',[],'YTick',[],
set(a,'YDir','normal')
a.XTick = 1:size(x,1);
a.YTick = 1:size(x,1);
xlabel('Model to be Fit');
ylabel('Model Data');
title(yl);