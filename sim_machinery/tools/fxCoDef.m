function CoDef = fxCoDef(y,yhat)
R = (y - yhat);    % Residuals
SR = R.^2; % Squared Residuals
SSR = sum(SR); % sum of squared residuals

E = (y - mean(y));    % Errors
SE = E.^2; % Squared Errors
SSE = sum(SE); % sum of squared errors

CoDef = 1-(SSR/SSE);




