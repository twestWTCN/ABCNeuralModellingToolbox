function Sigma = shrinkCorr(Sigma, alpha)
% SHRINKCORR  Shrink the correlation structure of Sigma toward the identity.
% Decomposes Sigma = D^{1/2} * Corr * D^{1/2}, shrinks off-diagonal
% elements of Corr by factor (1-alpha), then reconstructs.  Marginal
% variances (diagonal of Sigma) are unchanged.
stds  = sqrt(diag(Sigma));
Corr  = Sigma ./ (stds * stds');           % normalised correlation matrix
Corr  = (1 - alpha) * Corr + alpha * eye(size(Corr));  % shrink toward I
Sigma = Corr .* (stds * stds');            % restore original variances
Sigma = (Sigma + Sigma') / 2;             % enforce symmetry
