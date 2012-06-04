% laguerreBasisSetup.m
%
% Calculates basis set of weighted continuous Laguerre polynomials
% Note: returns oversampled signals
%
% optional: convolve Laguerre basis with input signal to produce signal Xi;
%       returned Xi is not oversampled.
%
% INPUTS:
%   p  = laguerre degree (p+1 basis functions = order 0:p);
%   a  = laguerre alpha (scaling);
%   len= basis length -> units == seconds;
%   Fs = sampling frequency
%   os = oversampling rate
%   is = input stimulus signal *** optional
%
% OUTPUTS:
%   weighted_laguerre_os = oversampled weighted continuous laguerre basis, 
%                           size = (len*os*Fs , p+1)
%   Xi = convolution of laguerre basis with input signal.
%
% 2012-01-12 Ben Cassidy


function [weighted_laguerre_os, Xi] = laguerreBasisSetup(p,a,len,Fs,os,is)

t = 0:1/(os*Fs):len;
t=t(:);

gen_laguerre = zeros(length(t),p+1);
for m=0:p
    % Using mfun is faster, but not all users have symbolic toolbox.
%     gen_laguerre(:,m+1) = mfun('L',m,a,t);
    for idx = 0:m
        gl_buf = ((-1)^idx)*nchoosek(m+a,m-idx)*(t.^idx)*(1/factorial(idx));
        gen_laguerre(:,m+1) = gen_laguerre(:,m+1) + gl_buf;
    end
end

weight = exp(-t).*t.^(a);
weightmatrix = repmat((weight), 1,p+1);
weighted_laguerre_os = (gen_laguerre.*weightmatrix);

% in case we didn't specify an input stimulus
if (~exist('is','var') || isempty(is)),Xi = [ ]; return, end 

% otherwise, convolve input stimulus with laguerre basis
is = is(:); % col vec
% is_osMat = repmat(is.',os,1);

is_osMat = is(:, ones(1,os)).';

is_os = is_osMat(:);
osNFFT = (2^nextpow2(max(numel(t),numel(is_os))));
is_k_os = (1/sqrt(osNFFT))*fft(is_os, osNFFT);
weighted_laguerre_k_os = (1/sqrt(osNFFT))*fft(weighted_laguerre_os,osNFFT);

Xi_k_os = weighted_laguerre_k_os.*repmat(is_k_os,1,p+1);
Xi_os = (sqrt(osNFFT))*ifft(Xi_k_os, 'symmetric');
Xi = Xi_os(1:os:end,:);
Xi = Xi(1:length(is),:);
end