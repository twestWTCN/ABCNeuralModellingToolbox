function R2 = fxSEucDist(y,yhat)
% DOI: 10.1214/16-BA1002
E = (y - yhat);    % Errors
E = E./std(yhat);
SE = sqrt(sum(E.^2)); 
R2 = -SE;

