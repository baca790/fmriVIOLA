% LM_NL.m
%
% Class for Lagrange Multiplier test -> fMRI applications
% Non-linear Hemodynamic Response model violation test.
%
% =========================================================================
% Use this to test for violation of the linaerity assumption for the
% hemodynamic response function, for each voxel time-series of an fMRI data
% volume. 
%
% USAGE: e.g. for test object "LM":
%
%   1)  Create a 'Non-Linearity test' object:
%   LM = LM_NL(data, input_stim_init, TR)
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


classdef LM_NL < LM_test_fmri
    
    properties
        len_HRF = 30;           % Length of HRF model (seconds);
        zeta_ab                 % Linear signal components
    end
    
    methods
        function lm_nl = LM_NL(data, input_stim_init, TR, extraDriftVars)
            if nargin < 4
                extraDriftVars = [];
                if nargin ~=3
                    data = [];
                    input_stim_init = [];
                    TR = [];
                end
            end
            lm_nl = lm_nl@LM_test_fmri(data, input_stim_init, TR, extraDriftVars);
        end
        
        function lm_nl = preLoopSetup(lm_nl)
            NFFT = lm_nl.NFFT;
            L = lm_nl.L;
            
            t = 0:lm_nl.TR:lm_nl.len_HRF;
            t = t(:);
            
            % Nonlinear model from (Purdon, et al., 2001)
            d_a = 1.5;
            d_b = 12;

            g_a_t = ((1 - exp(-1/d_a)).^2).*(t+1).*exp(-t./d_a);
            g_a_t = g_a_t./sum(g_a_t);
            g_a = timeFreqSig(g_a_t,'time', NFFT, L);
            g_b_t = (1 - exp(-1./d_b)).*exp(-t./d_b);
            g_b = timeFreqSig(g_b_t,'time', NFFT, L);            

            Xi.a = timeFreqSig(g_a.k.*lm_nl.input_stim.k, 'freq', NFFT, L);
            Xi.b = timeFreqSig(g_b.k.*lm_nl.input_stim.k,'freq', NFFT, L);
            Xi.c = timeFreqSig(Xi.b.t.*Xi.a.t, 'time', NFFT, L);

            % define zeta as the non-interaction signal components
            lm_nl.zeta_ab = timeFreqSig([Xi.a.t, Xi.b.t], 'time', NFFT, L);
            
            % assuming no temporal delay
            zeta_c_k = Xi.c.k;
            % calculate regressor
            lm_nl.z_k = zeta_c_k.';

            lm_nl = setupPolyDrift(lm_nl);
            
        end
        
        function lm_nl = inLoopSetup(lm_nl)
            
            BIC = NaN*zeros(lm_nl.ordermax_drift,...
                            lm_nl.ordermax_AR);
            for o_polydrift = 1:lm_nl.ordermax_drift
                drift.t = [lm_nl.polydrift.t(:,1:o_polydrift), ...
                    lm_nl.extraDriftVars.t];
                regMatx = [lm_nl.zeta_ab.t drift.t]; % just H0 since there 
                                                     % is no order
                                                     % selection for H1 in
                                                     % this case
                b = regMatx\lm_nl.y.t;
                v = lm_nl.y.t - regMatx*b;
                v_mc = v - mean(v);
                if ~isempty(lm_nl.ARpreallocated);
                    range_AR = lm_nl.ARpreallocated;
                else
                    range_AR = 0:lm_nl.ordermax_AR;
                end
                for o_AR = range_AR;
                    if o_AR ==0
                        nVar = ((lm_nl.y.t-regMatx*b).'*(lm_nl.y.t-regMatx*b))/(lm_nl.L - size(regMatx,2));
                    else
                        [~, nVar] = arburg(v_mc ,o_AR); % AR model
                    end
                    BIC(o_polydrift, o_AR+1) = ...
                        lm_nl.L*log(nVar) + log(lm_nl.L)*(o_polydrift+o_AR);
                end
            end
            % extract minimum BIC point here
            [~,BICminIdx] = min(BIC(:));
            
            % set orders based on minimum BIC
            [order_polydrift, lm_tv.order_AR] = ...
                ind2sub([lm_nl.ordermax_drift,...
                         lm_nl.ordermax_AR+1], BICminIdx);
            drift.t = [lm_nl.polydrift.t(:,1:order_polydrift), ...
                lm_nl.extraDriftVars.t];
            drift.k = [lm_nl.polydrift.k(:,1:order_polydrift), ...
                lm_nl.extraDriftVars.k];
                        
            lm_nl.X.t = [lm_nl.zeta_ab.t ,drift.t];
            lm_nl.X.k = [lm_nl.zeta_ab.k ,drift.k];
            
            lm_nl.x_k = [lm_nl.zeta_ab.k, drift.k].';
                        
%             if ~isempty(lm_nl.presetPolyDrift)
%                 order_polydrift = lm_nl.presetPolyDrift;
%             else
%                 order_polydrift = getRegressionModelOrder(lm_nl.y.t, lm_nl.polydrift.t);
%             end
%             drift.t = [lm_nl.polydrift.t(:,1:order_polydrift), ...
%                         lm_nl.extraDriftVars.t];
%             drift.k = [lm_nl.polydrift.k(:,1:order_polydrift), ...
%                         lm_nl.extraDriftVars.k];
%             lm_nl.DEBUG.order_drift = order_polydrift;
%             
%             lm_nl.X.t = [lm_nl.zeta_ab.t ,drift.t];
%             lm_nl.X.k = [lm_nl.zeta_ab.k ,drift.k];
%             
%             lm_nl.x_k = [lm_nl.zeta_ab.k, drift.k].';
        end
        
    end
end