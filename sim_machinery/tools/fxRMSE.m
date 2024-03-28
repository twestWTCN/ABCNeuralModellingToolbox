function RMSE = fxRMSE(y,yhat)

E = (y - yhat);    % Errors
SE = E.^2; % Squared Error
RMSE = sqrt(mean(SE));
RMSE = -RMSE;

% SSE = sum(SE); %negative sum of squared errors

