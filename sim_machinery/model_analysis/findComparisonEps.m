function [R P] = findComparisonEps(R,r2bank,exc)

prct = 50;
r2bankcat = horzcat(r2bank{:});
R.modcomp.modEvi.epspop = prctile(r2bankcat,prct); % threshold becomes median of model fits
for i = unique(exc)
    r2bankcat = horzcat(r2bank{exc==i});
    R.modcomp.modEvi.epspop(i) = prctile(r2bankcat,prct); % threshold becomes median of model fits
    for modID = find(exc==1)
        P(modID) = sum(r2bank{modID}>R.modcomp.modEvi.epspop(i))/size(r2bank{modID},2);
    end
end