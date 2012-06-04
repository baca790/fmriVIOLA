% LM_TV.m
%
% Class for Lagrange Multiplier test -> fMRI applications
% Time-Varying Hemodynamic Response model violation test.
%
% =========================================================================
% Use this to test for violation of the stationary-time assumption for the
% hemodynamic response function, for each voxel time-series of an fMRI data
% volume. 
%
% USAGE: e.g. for test object "LM":
%
%   1)  Create a 'time-varying HRF test' object:
%   LM = LM_TV(data, input_stim_init, TR)
%       data            = type=double  : dim (#TimePoints x #Voxels)
%                         -- fMRI data
%       input_stim_init = type=logical : dim (#TimePoints x 1)
%                         -- input stimulus signal
%       TR              = type=double  : scalar
%                         -- Sampling repetition time (seconds)
%
%   2)
%   LM = LM.apply_test_across_data;
%
%   3)
%   Test statistic results are stored in LM.Vt
%
%   Test statistic values in LM.Vt indicate model violation at
%   those voxels, above the threshold of the chi^2 distribution defined by
%   p_value and model dimensions.
%
% =========================================================================
%
%
%
%
% Ben Cassidy 2012-03-31

classdef LM_TV < LM_test_fmri
    
    properties
        ordermax_LBF = 13;       % laguerre basis regressor order
        ordermax_m = 5;          % max order of time-varying model
        a_laguerre = 5;          % Laguerre-basis scaling parameter
        len_HRF = 30;            % Length of HRF model (seconds);        
        Z_full;                  % Anomaly regression matrix
        Xi;                      % convolution of input stimulus with 
                                 % continuous Laguerre polynomial basis set
                                 % (L x ordermax_LBF).
        
        
    end
    
    methods
        function lm_tv = LM_TV(data, input_stim_init, TR, extraDriftVars)
            if nargin < 4
                extraDriftVars = [];
                if nargin ~= 3
                    data = [];
                    input_stim_init = [];
                    TR = [];
                end
            end
            lm_tv = lm_tv@LM_test_fmri(data, input_stim_init, TR, extraDriftVars);
        end
        
        function lm_tv = preLoopSetup(lm_tv)
            lm_tv = setupPolyDrift(lm_tv);
              
            [~, Xi_t] = laguerreBasisSetup(...
                    lm_tv.ordermax_LBF, ...
                    lm_tv.a_laguerre,...
                    lm_tv.len_HRF,...
                    1./lm_tv.TR,...
                    lm_tv.os_rate,...
                    lm_tv.input_stim.t);
            lm_tv.Xi = timeFreqSig(Xi_t, 'time', lm_tv.NFFT, lm_tv.L);
                                
            Z_full_t = fullTimeVaryingHRFvariables(lm_tv);
            lm_tv.Z_full = timeFreqSig(Z_full_t, 'time', lm_tv.NFFT, lm_tv.L);
        end
        
        function lm_tv = inLoopSetup(lm_tv)
            if ~isempty(lm_tv.presetPolyDrift)
                order_polydrift = lm_tv.presetPolyDrift;
            else
                order_polydrift = getRegressionModelOrder(lm_tv.y.t, lm_tv.polydrift.t, 1);
            end
            
            drift.t = [lm_tv.polydrift.t(:,1:order_polydrift), ...
                        lm_tv.extraDriftVars.t];
            drift.k = [lm_tv.polydrift.k(:,1:order_polydrift), ...
                        lm_tv.extraDriftVars.k];
            
            y_nodrift = lm_tv.y.t - drift.t*(drift.t\lm_tv.y.t);
            order_LBF = getRegressionModelOrder(y_nodrift, lm_tv.Xi.t);
%             lm_tv.DEBUG.order_LBF = order_LBF;
            
            order_TV = getRegressionModelOrder(y_nodrift, lm_tv.Z_full.t, lm_tv.ordermax_LBF+1, order_LBF);
%             lm_tv.DEBUG.order_TV = order_TV;
            
            lm_tv.X.t = [lm_tv.Xi.t(:,1:order_LBF) drift.t];
            lm_tv.X.k = [lm_tv.Xi.k(:,1:order_LBF) drift.k];
            lm_tv.x_k = [lm_tv.Xi.k(:,1:order_LBF) drift.k].';

            lm_tv.z_k = [];
            for index = 1:order_TV
            lm_tv.z_k = [lm_tv.z_k, lm_tv.Z_full.k(:,(1:order_LBF)+(index-1)*(lm_tv.ordermax_LBF+1))];
            end
            lm_tv.z_k = lm_tv.z_k.';
        end

        function Z_fullTVHRF = fullTimeVaryingHRFvariables(lm_tv)
            m = lm_tv.ordermax_m;
            L = lm_tv.L;
            time = (1:L)';
            
            w = 2*cos((1:m)*pi/(2*m+1));
            phi = zeros(L,m);
            for k=1:m
                phi(:,k) = (w(k)*cos(pi*k*time/L)).';
            end
            
            Xi_dim = size(lm_tv.Xi.t,2);
            Z = zeros(L, Xi_dim*m);
            for n=1:m
                Z(:,((n-1)*Xi_dim +1):((n-1)*Xi_dim)+Xi_dim) = repmat(phi(:,n), 1, Xi_dim).*(lm_tv.Xi.t);
            end
            
            Z_fullTVHRF = Z;
            
        end
        
        
        
    end
end
