% LM_DG.m
%
% Class for Lagrange Multiplier test -> fMRI applications
% Double Gamma (opt: + temporal derivative) model violation test
%
% =========================================================================
% Use this to test for violation of the Double Gamma hemodynamic response 
% function (optional: plus the Double Gamma temporal derivative as 
% specified in SPM or FSL software) for each voxel time-series of an fMRI
% data volume.
% 
% USAGE: e.g. for test object "LM":
%   
%   1)
%   LM = LM_DG(data, input_stim_init, TR, DG_DERIV);
%       data            = type=double  : dim (#TimePoints x #Voxels)
%                         -- fMRI data
%       input_stim_init = type=logical : dim (#TimePoints x 1)
%                         -- input stimulus signal
%       TR              = type=double  : scalar
%                         -- Sampling repetition time (seconds)
%       DG_DERIV        = type=string  : 'SPM', 'FSL', false
%                         -- optional, test validity of DG+derivative model
%   
%   2)
%   LM = LM.apply_test_across_data;
%
%
%   3)
%   Results are stored in LM.Vt
%
%   Test statistic values in LM.Vt indicate model violation at
%   those voxels, above the threshold of the chi^2 distribution defined by
%   p_value and model dimensions.
%
% =========================================================================
%
% Ben Cassidy 2012-03-29

classdef LM_DG < LM_test_fmri
    
    properties
        TEST_DERIVATIVE;         % test DG derivative:'SPM'||'FSL'||'false'
        ordermax_LBF = 14;       % laguerre basis regressor order
        Xi                       % convolution of input stimulus with 
                                 % continuous Laguerre polynomial basis set
                                 % (L x ordermax_LBF).
        rho                      % convolution of Double-Gamma (and 
                                 % derivative function if specified)  with
                                 % input stimlulus signal.
        a_laguerre = 5;          % Laguerre-basis scaling parameter for 
                                 % standard Double-Gamma test. 
                                 % This changes if testing the 'Double 
                                 % Gamma plus derivative' model; uses the
                                 % highest valid order for each test to
                                 % preserve degrees of freedom.
        len_HRF = 30;            % Length of HRF model (seconds);
    end
    
    methods
        % class constructor
        function lm_dg = LM_DG(data, input_stim_init, TR, DG_DERIV, extraDriftVars)
            if nargin < 5
                if nargin < 4
                    if nargin ~= 3
                        data = [];
                        input_stim_init = [];
                        TR = [];
                    end
                    DG_DERIV = false;
                end
                extraDriftVars = [];
            end
            lm_dg = lm_dg@LM_test_fmri(data, input_stim_init, TR, extraDriftVars);
            assert(isequal(DG_DERIV,false) || ...
                strcmp(DG_DERIV,'SPM') || ...
                strcmp(DG_DERIV,'FSL'));
            lm_dg.TEST_DERIVATIVE = DG_DERIV;
        end

        % preLoopSetup initialises most signals
        function lm_dg = preLoopSetup(lm_dg)
            lm_dg = setupPolyDrift(lm_dg);
            lm_dg = doubleGammaDerivSetup(lm_dg); % Init rho
            [~, Xi_t] = laguerreBasisSetup(...
                    lm_dg.ordermax_LBF, ...
                    lm_dg.a_laguerre,...
                    lm_dg.len_HRF,...
                    1./lm_dg.TR,...
                    lm_dg.os_rate,...
                    lm_dg.input_stim.t);
            lm_dg.Xi = timeFreqSig(Xi_t, 'time', lm_dg.NFFT, lm_dg.L);
        end
        
        % Signals which require updating for every voxel, i.e. model order
        % selection.
        function lm_dg = inLoopSetup(lm_dg)
            
            BIC = NaN*zeros(lm_dg.ordermax_drift,...
                            lm_dg.ordermax_LBF+1,...
                            lm_dg.ordermax_AR+1);
            for o_polydrift = 1:lm_dg.ordermax_drift
                drift.t = [lm_dg.polydrift.t(:,1:o_polydrift), ...
                    lm_dg.extraDriftVars.t];
                for o_LBF = 1:lm_dg.ordermax_LBF+1
                    LBF.t = lm_dg.Xi.t(:,1:o_LBF);
                    
                    regMatx = [drift.t LBF.t];
                    b = regMatx\lm_dg.y.t;
                    v = lm_dg.y.t - regMatx*b;
                    v_mc = v - mean(v);
                    
                    if ~isempty(lm_dg.ARpreallocated);
                        range_AR = lm_dg.ARpreallocated;
                    else
                        range_AR = 0:lm_dg.ordermax_AR;
                    end
                    for o_AR = range_AR;
                        if o_AR ==0
                            nVar = ((lm_dg.y.t-regMatx*b).'*(lm_dg.y.t-regMatx*b))/(lm_dg.L - size(regMatx,2));
                        else
                            [~, nVar] = arburg(v_mc ,o_AR); % AR model
                        end
                        BIC(o_polydrift, o_LBF, o_AR+1) = ...
                            lm_dg.L*log(nVar) + log(lm_dg.L)*(o_polydrift+o_LBF+o_AR);
                    end
                end
            end
            
            % extract minimum BIC point here
            [~,BICminIdx] = min(BIC(:));
            
            % set orders based on minimum BIC
            [order_polydrift, order_LBF, lm_tv.order_AR] = ...
                ind2sub([lm_dg.ordermax_drift,...
                         lm_dg.ordermax_LBF+1,...
                         lm_dg.ordermax_AR+1], BICminIdx);
            drift.t = [lm_dg.polydrift.t(:,1:order_polydrift), ...
                lm_dg.extraDriftVars.t];
            drift.k = [lm_dg.polydrift.k(:,1:order_polydrift), ...
                lm_dg.extraDriftVars.k];
            
            lm_dg.X.t = [lm_dg.rho.t, drift.t];
            lm_dg.X.k = [lm_dg.rho.k, drift.k];
            lm_dg.x_k = [lm_dg.rho.k, drift.k].';
            lm_dg.z_k = lm_dg.Xi.k(:,1:order_LBF).';
            
%             if ~isempty(lm_dg.presetPolyDrift)
%                 order_polydrift = lm_dg.presetPolyDrift;
%             else
%                 order_polydrift = getRegressionModelOrder(lm_dg.y.t, lm_dg.polydrift.t, 1);
%             end
%             
%             drift.t = [lm_dg.polydrift.t(:,1:order_polydrift), ...
%                         lm_dg.extraDriftVars.t];
%             drift.k = [lm_dg.polydrift.k(:,1:order_polydrift), ...
%                         lm_dg.extraDriftVars.k];
%           
%             y_nodrift = lm_dg.y.t - drift.t*(drift.t\lm_dg.y.t);
%             order_LBF = getRegressionModelOrder(y_nodrift, lm_dg.Xi.t);
% %             lm_dg.DEBUG.order_LBF = order_LBF;
% %             lm_dg.DEBUG.order_drift = order_polydrift;
%             
%             lm_dg.X.t = [lm_dg.rho.t, drift.t];
%             lm_dg.X.k = [lm_dg.rho.k, drift.k];
%             lm_dg.x_k = [lm_dg.rho.k, drift.k].';
%             lm_dg.z_k = lm_dg.Xi.k(:,1:order_LBF).';
        end
        
        % setup convolved HRF signal
        function lm_dg = doubleGammaDerivSetup(lm_dg)
            
            t = 0:lm_dg.TR:lm_dg.len_HRF;
            t = t(:);
            
            % standard SPM, FSL double gamma fn
            dg_HRF_t= ((t.^5.*exp(-t))/gamma(6) ...
                - (1.0/6)*(t.^15.*exp(-t))/gamma(16));
            dg_HRF_norm_t = dg_HRF_t./sum(dg_HRF_t);
            dg_HRF = timeFreqSig(dg_HRF_norm_t, 'time', lm_dg.NFFT, lm_dg.L);
            rho_k = dg_HRF.k .* lm_dg.input_stim.k;
            
            % temporal derivative - (Cassidy & Solo, 2012, ISBI)
            switch upper(lm_dg.TEST_DERIVATIVE)
                case 'FSL'
                    dg_derivative_t = ((t.^4).*exp(-t)).*((5-t)./gamma(6) - ...
                        (15.*(t.^10) - t.^11)./(6*gamma(16)));
                    lm_dg.a_laguerre = 4;
                case 'SPM'
                    sum1=zeros(size(t));
                    sum2=zeros(size(t));
                    for k=0:5
                        sum1 = sum1+ nchoosek(5,k).*((-1).^(5-k)).*(t.^k).*exp(1);
                    end
                    sum1(t<1)=0;
                    for k=0:15
                        sum2 = sum2+ nchoosek(15,k).*((-1).^(15-k)).*(t.^k).*exp(1);
                    end
                    sum2(t<1)=0;
                    
                    dg_derivative_t = exp(-t).*((((t.^5) - (sum1))./(gamma(6))) ...
                        - (((t.^15) - (sum2))./(6*gamma(16))));
                    lm_dg.a_laguerre = 0;
                otherwise
                    %
                    lm_dg.a_laguerre = 5;
            end
            
            try
                dg_derivative_norm_t = dg_derivative_t./sum(dg_derivative_t);
                dg_derivative = timeFreqSig(dg_derivative_norm_t, 'time', lm_dg.NFFT, lm_dg.L);
                rho_k = [rho_k, (dg_derivative.k.*lm_dg.input_stim.k)];
            catch
            end
            
            lm_dg.rho = timeFreqSig(rho_k, 'freq', lm_dg.NFFT, lm_dg.L);
        end
    end

end