% bonferroni_chi2.m
%
% Calculates threshold from chi-square distribution with Bonferroni 
% correction.
%
% input is either:
%   (x,p,p_value), x = number of comparisons, 
%                  p = degrees of freedom,
%                  p_value = probability
% or
%   (x,y,z, p, p_value), [x,y,z] is image volume dimensions if you are lazy
%                  but this is probably not useful for
%                  fmri since volume usually includes non-brain voxels
%                  which were not tested.
%
% Ben Cassidy 2009-01-30
% changed 29/10/10 to allow any % c.i.

function threshold = bonferroni_chi2(varargin)

if (nargin == 3)
    x = varargin{1};
    p = varargin{2};
    pvalue = varargin{3};
elseif (nargin == 5)
    x = varargin{1}*varargin{2}*varargin{3};
    p = varargin{4};
    pvalue = varargin{5};
else
    error('Number of input arguments must be 2  or 4, type "help bonferroni_95pc_chi2" for details');
end    

threshold =  chi2inv(1-(pvalue/x), p);