% offset_drift_orthogonal.m
%
% Construct polynomial drift regressor matrix, shape =(length x (order +1))
%
% Optional input 'noConst' returns basis without constant (dc) vector
%
% Ben Cassidy 2012-03-29


function [basis] =  offset_drift_orthogonal(order, length, noConst)

T=length;

if (mod(T,2) ==1)
    t0 = ones(T,1);
    t1 = (((1:T)-T/2).');
    t2 = (((1:T)-T/2).').^2;
    t3 = (((1:T)-T/2).').^3;
else
    t0 = ones(T,1);
    t1 = (((1:T)-(T-1)/2).');
    t2 = (((1:T)-(T-1)/2).').^2;
    t3 = (((1:T)-(T-1)/2).').^3;
end

% make b0 independent of first order term
t1_res = t1 - t0;

% estimate 2nd order error
y=t2;
X = [t0 t1];
alpha= X\y;
t2_reg = X*alpha;

t2_res = t2 - t2_reg;

% repeat for 3rd order
y=t3;
X = [t0 t1_res t2_res];
alpha = X\y;
t3_reg = X*alpha;

t3_res = t3 - t3_reg;

switch order
    case 0
        basis = [t0];
    case 1
        basis = [t0 t1_res];
    case 2
        basis = [t0 t1_res t2_res];
    case 3
        basis = [t0 t1_res t2_res t3_res];
    otherwise
        %default use O(3)
        basis = [t0 t1_res t2_res t3_res];
        printf('Using default o(3) drift basis set')
end

basis = basis./repmat(max(basis), T,1);
try basis(:,2:2:end) = -basis(:,2:2:end); end 
try if noConst, basis = basis(:,2:end); end, end