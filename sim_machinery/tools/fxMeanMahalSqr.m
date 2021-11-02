function R2 = fxMeanMahalSqr(y,yhat)
D = mahal(y,yhat);
R2 = mean(D.^2);
R2 = -R2; % negate, as algorithm is maximizing
