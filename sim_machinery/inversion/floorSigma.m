function Sigma = floorSigma(Sigma, priorSigma, alpha)
% FLOORSIGMA  Apply a hard eigenvalue floor to a covariance matrix.
%
%   Sigma = floorSigma(Sigma, priorSigma, alpha)
%
%   Sets a lower bound on each eigenvalue of Sigma equal to
%   alpha * (smallest eigenvalue of priorSigma).  This ensures the
%   posterior Gaussian approximation never becomes more concentrated than
%   alpha times the sharpest prior direction, without pulling the estimate
%   back toward the prior (as additive mixing does).
%
%   alpha = 0.01 (default) allows 100x concentration relative to prior.
%   Increase alpha to slow posterior focusing; decrease to allow more.

priorEigs   = eig(priorSigma);
eigFloor    = alpha * min(priorEigs(priorEigs > 0));

[V, D]   = eig(Sigma);
d        = diag(D);
d        = max(d, eigFloor);
Sigma    = V * diag(d) * V';
Sigma    = (Sigma + Sigma') / 2; % enforce symmetry after floating-point ops