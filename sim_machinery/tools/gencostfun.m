function [R2] = gencostfun(y,yhat)
E = (y - yhat);    % Errors
SE = E.^2; %(y - yhat).^2   % Squared Error

Epr = (y - mean(y));
SEpr = Epr.^2;
MSE = sum(SE)/sum(SEpr); %mean((y - yhat).^2)
RMSE = sqrt(MSE); %sqrt(mean((y - yhat).^2));  
