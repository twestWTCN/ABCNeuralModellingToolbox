function fitness = sampleFitness(indices,R,parBank,pIndMap,pOrg)
selected = find(indices);
% Ensure the minimum sample size is met
if numel(selected) < R.SimAn.minRank
    % Penalize configurations with fewer than minSampleSize samples
    fitness = inf;
else
    % Calculate mean error of selected samples
    parOptBank = parBank(:, selected);
    meanError = mean(parOptBank(end,:));
    % Calculate KL divergence of the selected parameter set
    Mfit.prior = R.Mfit.prior;
    

    s = parOptBank(end,:);
    xs = parOptBank(pIndMap,:);
    W = ((s(end,:)-1).^-1);
    W = W./sum(W);
    Ws = repmat(W,size(xs,1),1); % added 03/2020 as below wasnt right dim (!)
    Mfit.Mu = wmean(xs,Ws,2);
    Mfit.Sigma = weightedcov(xs',W);

        [~,DKL,R] = KLDiv(R,Mfit,pOrg,[],0);

        % Objective function to minimize (negative)
        fitness = -(R.SimAn.scoreweight(1)*meanError - R.SimAn.scoreweight(2)*DKL);

end
% figure(100)
% imagesc(indices)
%  title(sprintf('%.3f',fitness))
%  drawnow
