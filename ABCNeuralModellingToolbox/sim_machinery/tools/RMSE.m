function [MSE RMSE] = RMSE(y,yhat)
E = (y - yhat);    % Errors
SE = E.^2; %(y - yhat).^2   % Squared Error
MSE = mean(SE); %mean((y - yhat).^2)   % Mean Squared Error
RMSE = sqrt(MSE); %sqrt(mean((y - yhat).^2));  % Root Mean Squared Error
