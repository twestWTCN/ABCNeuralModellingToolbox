function SSE = RMSE_scaled(y,yhat)
E = (y - yhat);    % Errors
SE = E.^2; %(y - yhat).^2   % Squared Error
SSE = sum(SE); %/sum(SEpr); %mean((y - yhat).^2)

