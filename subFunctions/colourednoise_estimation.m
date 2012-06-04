% colourednoise_estimation.m
%
% [beta_est, err_sum] = colourednoise_estimation(y, X, b)
%
% Estimates new parameters for AR(n) + WN noise signal
% inputs (y, X, b) for regression model y=X*b + v
%
% Uses AR(p) model for p= order calculated using BIC
% unless (optional) p=exactARorder is specified
%
% Operation:
% 1) uses input beta to calculate error
% 2) fits AR(p) to error to get AR parameter 'phi'
% 3) filters y by AR parameter to get y_est
% 4)    "    X    "   "    "       "  X_est
% 5) calculate new value of beta by least squares regression
%
% Ben Cassidy - 2009/03/04

function [b_out, order_out, BIC, y_out, X_out] = colourednoise_estimation(y, X, b, exactARorder)
try
    if exactARorder ==0
        b_out = b;
        order_out = 0;
        BIC = [];
        y_out = y;
        X_out = X;
        return
    end
catch % exactARorder is not specified
end

n=length(y);

maxOrder = 7;
order =0;
order_out=0;
BIC = NaN*zeros(maxOrder+1, 1);

residual_var = ((y-X*b).'*(y-X*b))/(n - size(X,2));

BIC(1) = n*log(residual_var); % zero-order
BICmin = BIC(1);

b_out = b;
y_out = y;
X_out = X;

while(order < maxOrder) % check all orders, it might have local minima
    order = order+1;
    
    % AR(p) estimation
    % ------------------------
    v = y - X*b;      % residual
    v_mc = v-mean(v); % mean-corrected
    
    [phi, nVar] = arburg(v_mc ,order); % AR model

    y_est = filtfilt(phi(:),1,y);  % signal estimated
    X_est = filtfilt(phi(:),1,X);  % parameter estimated
    b_est = X_est\y_est;
    % ------------------------
    
    % BIC order selection
    BIC(order+1) = n*log(nVar) + order*log(n);
    
    try
        % iff we only want the AR(p) estimate for prespecified p
        if (exactARorder == order)
            b_out = b_est;
            order_out = order;
            y_out = y_est;
            X_out = X_est;
            return,
        end
    catch
    end
    % else we just iterate through to find the ideal AR(p) order wrt BIC
    if BIC(order+1) < BICmin
        BICmin = BIC(order+1);
        b_out = b_est;
        order_out = order;
        y_out = y_est;
        X_out = X_est;
    end
end
end