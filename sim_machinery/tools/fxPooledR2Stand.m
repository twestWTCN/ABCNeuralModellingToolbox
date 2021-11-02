function R2 = fxPooledR2Stand(y,yhat)
E = (y - yhat);    % Errors
SE = sum(E.^2); % Squared Error
T = (y-mean(y)); 
ST = sum(T.^2);

R2 = 1-(SE/ST);

