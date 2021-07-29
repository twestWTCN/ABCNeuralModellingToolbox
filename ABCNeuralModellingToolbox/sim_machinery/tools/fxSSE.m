function SSE = fxSSE(y,yhat)
E = (y - yhat);    % Errors
SE = E.^2; % Squared Error
SSE = -sum(SE); %negative sum of squared errors

