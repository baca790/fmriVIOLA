% LM_test_fmri.m
%
% superclass for Lagrange Multiplier test -> fMRI applications
% superclass: LM_test
% subclasses: individual test definitions
%
% Do not call this class definition directly, instead use via individual
% test definitions
%
% 2012-03-29 Ben Cassidy

classdef LM_test_fmri < LM_test
    
    properties
        os_rate = 20;        % for generating accurate orthogonal basis expansions
        NFFT                 % number of FFT points
        V                    % Test statistic - all voxels
        Vt                   % Test statistic (threshold limited)
        numVoxels            %
        data                 % temporal mean-corrected fMRI data, matrix size (L x numVoxels)
        datamean             % temporal mean of raw data
        input_stim_init      % BOOL: column vector, length == L
        input_stim           % struct: .t = time signal, .k = frequency sig 
        TR                   % fMRI sampling time (seconds)
        p_value = 0.01       %
        threshold            % calculated from chi^2 distribution
        numComparisons = 1;  % initialise - adjust in subclasses if needed
        polydrift            % drift regressors, (L x ordermax_drift)
                             % struct: .t = time signal, .k = freq sig
        ordermax_drift = 3;  % Default maximum polynomial drift order
        ordermax_AR = 3;     % Default maximum autoregressive model order
        order_AR             % ?private: AR model order in use
        presetPolyDrift =[]; % integer: default is to ignore this and 
                             % choose polynomial drift order using BIC for 
                             % an accurate drift model fit.
                             % Otherwise use this to specify integer 
                             % from 0 (just dc) to 3 (all up to 3rd 
                             % order polynomial drift)
        extraDriftVars       % Default: disabled, otherwise specify this as
                             % other drift regressors to include (e.g.
                             % movement regressors or a different drift
                             % basis.
                             % At the moment, this disables the automatic
                             % polynomial drift order selection, defaults
                             % to the maximum (i.e. ordermax_drift) and 
                             % includes all polynomial drift regressors 
                             % unless the presetPolyDrift variable is set
                             % to a different order.
        ARpreallocated = []; % integer: Autoregressive model order 
                             % (default empty: choose automatically), 
        DEBUG
    end

    methods
        % class constructor
        function lm_fmri = LM_test_fmri(data, input_stim_init, TR, extraDriftVars)
            if nargin < 4
                extraDriftVars =[];
                if nargin ~= 3
                    input_stim_init = [];
                    data = [];
                    TR = [];
                end
            end
            
            lm_fmri = lm_fmri@LM_test( ); % placeholder for other init.
            lm_fmri.data = LM_test_fmri.dataReshape(double(data));
            [lm_fmri.L, lm_fmri.numVoxels] = size(lm_fmri.data);
            lm_fmri.NFFT = (2^nextpow2(lm_fmri.L));
            lm_fmri.input_stim = timeFreqSig(+input_stim_init, 'time', ...
                                            lm_fmri.NFFT, lm_fmri.L);
            lm_fmri.TR = TR;
            lm_fmri.V = zeros(lm_fmri.numVoxels,1);
            lm_fmri.Vt = zeros(lm_fmri.numVoxels,1);
            lm_fmri.extraDriftVars = timeFreqSig(+extraDriftVars, ...
                                          'time', lm_fmri.NFFT, lm_fmri.L);
        end
        
        % Call this function to calculate the test statistic across all
        % voxels
        function lm_fmri = apply_test_across_data(lm_fmri)
            % preLoopSetup should contain most model setup calculations 
            lm_fmri = preLoopSetup(lm_fmri);
            L = lm_fmri.L;
            NFFT_buf = lm_fmri.NFFT;
            
            nonzeroVox = sum(abs(bsxfun(@minus, lm_fmri.data, ...
                                                mean(lm_fmri.data)))) > 0;
            
            % loop through voxels
            for n=1:lm_fmri.numVoxels
                y_t = lm_fmri.data(:,n);
                
                % do not process voxels containing null data
                if nonzeroVox(n)
                    lm_fmri.y = timeFreqSig(y_t, 'time', NFFT_buf, L);
                    
                    % inLoopSetup should only contain data-dependent
                    % calculations, e.g. model order selection. Preallocate
                    % outside this loop.
                    lm_fmri = inLoopSetup(lm_fmri);
                    
                    % fit model under null hypothesis, pre-AR
                    b =lm_fmri.X.t\lm_fmri.y.t;
                    
                    if lm_fmri.order_AR == 0
                        y_ar = lm_fmri.y.t;
                        X_ar = lm_fmri.X.t;
                        lm_fmri.b = b;
                    else
                        v = lm_fmri.y.t - lm_fmri.X.t*b; % residual
                        v_mc = v - mean(v);
                        
                        [phi, ~] = arburg(v_mc ,lm_fmri.order_AR); % AR model
                        y_ar = dcfilt(phi(:),1,lm_fmri.y.t);  % signal estimated
                        X_ar = dcfilt(phi(:),1,lm_fmri.X.t);  % regressor estimated
                        lm_fmri.b = X_ar\y_ar;
                    end

                                             
                    % AR corrected signals
                    lm_fmri.y = timeFreqSig(y_ar, 'time', NFFT_buf, L);
                    lm_fmri.X = timeFreqSig(X_ar, 'time', NFFT_buf, L);
                    
                    % Test statistic at voxel (n)
                    lm_fmri.V(n) = lagrange_multiplier_test_univariate(lm_fmri);
                end

                if (mod(n,500) ==0), disp([num2str(100*n/lm_fmri.numVoxels) ' %']), end
            end
            
            lm_fmri = apply_threshold_across_data(lm_fmri);
        end

        function lm_fmri = apply_threshold_across_data(lm_fmri)
            lm_fmri.numComparisons = length(find(lm_fmri.V));
            lm_fmri.threshold = bonferroni_chi2(...
                    lm_fmri.numComparisons, ...
                    size(lm_fmri.X.t,2), ...
                    lm_fmri.p_value);
            lm_fmri.Vt = lm_fmri.V;
            lm_fmri.Vt(lm_fmri.V < lm_fmri.threshold) = NaN;
        end
        
        function lm_fmri = setupPolyDrift(lm_fmri)
            polydrift_t = offset_drift_orthogonal(lm_fmri.ordermax_drift, lm_fmri.L);
            lm_fmri.polydrift = timeFreqSig(polydrift_t, 'time', lm_fmri.NFFT, lm_fmri.L);            
        end
        
    end
    
    methods (Static)
        function data = dataReshape(data)
            % Use this function to reshape / validate fmri data input
            % currently assumes data shape is (numTimeSamples x numVoxels)
        end
        
    end
end

function out = dcfilt(B, A, X)
% filtering w/ dc removal; this is a small hack to avoid having to deal with
% initial AR conditions, e.g. via backwards forecasting initial estimate
dc = (X(1,:));
X = bsxfun(@minus, X, dc);
out = filter(B, A, X);
out = bsxfun(@plus, out, dc);
end