% getRegressionModelOrder.m
%
%
% groupSize is an optional input; if regressorVariables is divided into groups (e.g. for
% coordinates or basis sets etc) it specifies the number of variable
% sub-components in each group.
% 
% If groupSize is specified, it returns [order] as the number of selected
% groups instead of the total number of variables 
% i.e. size(regressorVariables,2) = groupSize * (number of groups)
%
%
% Ben Cassidy 2012-03-29

function [modelOrder] = getRegressionModelOrder(data, regressorVariables, groupSize, subGroupPickElts)
n=length(data);
if (~exist('groupSize', 'var') || groupSize ==1 || ~exist('subGroupPickElts', 'var'))
    groupSize =1;
    subGroupPickElts = 1;
end

if groupSize ==1
    numberOfGroups =size(regressorVariables,2);
else
    numberOfGroups =size(regressorVariables,2)/groupSize;
end

BIC_order = zeros(numberOfGroups,1)*NaN;
for k=1:numberOfGroups

    if groupSize == 1
        X = regressorVariables(:,1:k);
    else
        X = []; %init.
        % lazy but it works and is understandable: iteratively select  and
        % concatenate the first columns of each group, then 2nd, 3rd etc.
        for index = 1:k
            X = [X, regressorVariables(:,(1:subGroupPickElts)+(index-1)*groupSize)];
        end
    end

%     X = regressorVariables(:,(1:k*groupSize));

    b = X\data;
    res = data- X*b;
    residual_var = (norm(res))^2/(n - k*subGroupPickElts);
    % change to resv = 1/(xxxx); resv = norm()^2*resv ..... speedup?
    
    BIC_order(k) = n*(log(residual_var)) + (k)*log(n);
end
    
[~,modelOrder] = min(BIC_order);
end