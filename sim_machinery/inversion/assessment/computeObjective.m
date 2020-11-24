function [ACC R2 DKL] = computeObjective(R,r2)
if ~isfield(R.SimAn,'scoreweight')
    R.SimAn.scoreweight = [1 0];
    warning('Combined score weightings not available, setting to default goodness of fit only')
end
R2 = R.SimAn.scoreweight(1)*(r2);
ACC = (R.SimAn.scoreweight(1)*(r2)) - (R.SimAn.scoreweight(2)*R.Mfit.DKL);
DKL = R.Mfit.DKL;