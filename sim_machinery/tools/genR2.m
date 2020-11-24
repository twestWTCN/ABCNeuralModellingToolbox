function [c] = genR2(y,yhat)
E = (y - yhat);    % Errors
SE = E.^2; %(y - yhat).^2   % Squared Error

Epr = (y - mean(y));
SEpr = Epr.^2;
c = sum(SE)/sum(SEpr); %mean((y - yhat).^2)

