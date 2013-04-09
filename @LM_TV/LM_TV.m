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
        ordermax_LBF = 12;       % laguerre basis regressor order
        ordermax_m = 4;          % max order of time-varying model
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
            
            % assuming no preset poly drift order
            BIC = NaN*zeros(lm_tv.ordermax_drift,...
                            lm_tv.ordermax_LBF+1,...
                            lm_tv.ordermax_m,...
                            lm_tv.ordermax_AR+1);
            for o_poly = 1:lm_tv.ordermax_drift
                drift.t = [lm_tv.polydrift.t(:,1:o_poly), ...
                                lm_tv.extraDriftVars.t];
                for o_LBF = 1:lm_tv.ordermax_LBF+1
                    LBF.t = lm_tv.Xi.t(:,1:o_LBF);
                    for o_TV = 1:lm_tv.ordermax_m
                        TV_cmpt.t = [];
                        for index = 1:o_TV
                            TV_cmpt.t = [TV_cmpt.t lm_tv.Z_full.t(:,(1:o_LBF)+(index-1)*lm_tv.ordermax_LBF+1)];
                        end
                        
                        %%% fit model with specified orders
                        %regression matrix for calculating model fit
                        %(unfortunately we have to do this under H_1)
                        regMatx = [drift.t LBF.t TV_cmpt.t];
                        
                        % intervening.....
                        % slight hack check for overfit which would lead to
                        % conditioning problems, collinear model etc.
                        % Shouldn't use these cases in valid analysis anyway.
                        if size(regMatx,2) >  (lm_tv.L / 10) % magic number, oh well...
                            break % BIC will be ignored for this combo.
                        end
                        
                        b = regMatx\lm_tv.y.t;
                        v = lm_tv.y.t - regMatx*b; % residual
                        v_mc = v- mean(v);
                        
                        
                        if ~isempty(lm_tv.ARpreallocated);
                            range_AR = lm_tv.ARpreallocated;
                        else
                            range_AR = 0:lm_tv.ordermax_AR;
                        end
                        for o_AR = range_AR;
                        
                            
                            % calculate BIC from model
                            if o_AR ==0
                                nVar = ((lm_tv.y.t-regMatx*b).'*(lm_tv.y.t-regMatx*b))/(lm_tv.L - size(regMatx,2));
                            else
                                [~, nVar] = arburg(v_mc ,o_AR); % AR model
                            end
                            BIC(o_poly, o_LBF, o_TV, o_AR+1) = ...
                                lm_tv.L*log(nVar) + log(lm_tv.L)*(o_poly+o_LBF+o_TV+o_AR);
                        end
                    end
                end
            end
            
            % extract minimum BIC point here
            [~,BICminIdx] = min(BIC(:));
            
            % set orders based on minimum BIC
            [order_polydrift, order_LBF, order_TV, lm_tv.order_AR] = ...
                ind2sub([lm_tv.ordermax_drift,...
                            lm_tv.ordermax_LBF+1,...
                            lm_tv.ordermax_m,...
                            lm_tv.ordermax_AR+1], BICminIdx);
            drift.t = [lm_tv.polydrift.t(:,1:order_polydrift), ...
                lm_tv.extraDriftVars.t];
            drift.k = [lm_tv.polydrift.k(:,1:order_polydrift), ...
                lm_tv.extraDriftVars.k];
            
            % set regressors based on min-orders
            lm_tv.X.t = [lm_tv.Xi.t(:,1:order_LBF) drift.t];
            lm_tv.X.k = [lm_tv.Xi.k(:,1:order_LBF) drift.k];
            lm_tv.x_k = [lm_tv.Xi.k(:,1:order_LBF) drift.k].';
            
            lm_tv.z_k = [];
            for index = 1:order_TV
                lm_tv.z_k = [lm_tv.z_k, lm_tv.Z_full.k(:,(1:order_LBF)+(index-1)*(lm_tv.ordermax_LBF+1))];
            end
            lm_tv.z_k = lm_tv.z_k.';
            
            
            
%             if ~isempty(lm_tv.presetPolyDrift)
%                 order_polydrift = lm_tv.presetPolyDrift;
%             else
%                 order_polydrift = getRegressionModelOrder(lm_tv.y.t, lm_tv.polydrift.t, 1);
%             end
%             
%             drift.t = [lm_tv.polydrift.t(:,1:order_polydrift), ...
%                         lm_tv.extraDriftVars.t];
%             drift.k = [lm_tv.polydrift.k(:,1:order_polydrift), ...
%                         lm_tv.extraDriftVars.k];
%             
%             y_nodrift = lm_tv.y.t - drift.t*(drift.t\lm_tv.y.t);
%             order_LBF = getRegressionModelOrder(y_nodrift, lm_tv.Xi.t);
% %             lm_tv.DEBUG.order_LBF = order_LBF;
%             
%             order_TV = getRegressionModelOrder(y_nodrift, lm_tv.Z_full.t, lm_tv.ordermax_LBF+1, order_LBF);
% %             lm_tv.DEBUG.order_TV = order_TV;
%             
%             lm_tv.X.t = [lm_tv.Xi.t(:,1:order_LBF) drift.t];
%             lm_tv.X.k = [lm_tv.Xi.k(:,1:order_LBF) drift.k];
%             lm_tv.x_k = [lm_tv.Xi.k(:,1:order_LBF) drift.k].';
% 
%             lm_tv.z_k = [];
%             for index = 1:order_TV
%             lm_tv.z_k = [lm_tv.z_k, lm_tv.Z_full.k(:,(1:order_LBF)+(index-1)*(lm_tv.ordermax_LBF+1))];
%             end
%             lm_tv.z_k = lm_tv.z_k.';
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
