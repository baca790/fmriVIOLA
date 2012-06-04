% timeFreqSig.m
%
% Input column-vector signal matrix
% Outputs a struct with the time-domain and frequency-domain signals of the
% input.  Outputs normalised fft and ifft values, assumes same for inputs.
%
%   IN
%
% sig       : double : time or frequency signal to convert
% type      : string : type of sig, either 'time' or 'freq'
% NFFT      : int : number of frequency samples
% L         : int : number of time samples
%
%   OUT
%
% outSig.t : double : time-domain signal
% outSig.k : double : frequency-domain signal
%
% Ben Cassidy 2012-03-29

function [outSig] = timeFreqSig(sig, type, NFFT, L)
switch lower(type)
    case 'time'
        outSig.t = sig;
        outSig.k = (1/sqrt(NFFT))*fft(sig, NFFT);
    case 'freq'
        outSig.t = (sqrt(NFFT))*ifft(sig, 'symmetric');
        outSig.t = outSig.t(1:L,:);
        outSig.k = sig;
    otherwise
        error('Invalid input type, should be "time" or "freq"')
end
end