function RMSE = fxzRMSE(y,yhat)
XStar = mean([y; yhat]);
SStar = std([y; yhat]);

y = (y-XStar)./SStar;
yhat = (yhat-XStar)./SStar;

E = (y - yhat);    % Errors
SE = E.^2; % Squared Error
RMSE = sqrt(mean(SE));
RMSE = -RMSE;

% SSE = sum(SE); %negative sum of squared errors

