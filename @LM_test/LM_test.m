% LM_test.m
%
% Superclass for Lagrange Multiplier test
%
% Do not call this class definition directly, instead use via individual
% test definitons -> this function is all tip and no iceberg.
%
% 2012-03-29 Ben Cassidy

classdef LM_test
    % General class for Lagrange Multiplier diagnostic tests.
    % Make subclasses for individual tests
    
    properties
        L                   % length of data signal
        V_univariate        % univariate test statistic
        y                   % data signal
        X                   % regression matrix: null hypothesis parameters
        z_k                 % alt regression matrix evaluated under H0
        x_k                 % null regression matrix evaluated under H0
        b                   % estimated parameters
    end

    methods
        % class constructor
        function lm = LM_test()
            % initialise everything in subclasses
        end
        
        % Test statistic calculation - univariate. 
        function V_univariate = lagrange_multiplier_test_univariate(lm)
            F_var = ((lm.y.t - lm.X.t*lm.b).'*(lm.y.t - lm.X.t*lm.b))/(lm.L-size(lm.X.t,2));
            
            % residuals; freq domain 'b' parameters are same as time domain
            v_k = (lm.y.k - lm.X.k*lm.b).';

            D_1 = real(lm.z_k*lm.z_k')/F_var;
            D_2 = real(lm.z_k*lm.x_k')/F_var;
            D_3 = real(lm.x_k*lm.x_k')/F_var;
            D_4 = real(lm.x_k*lm.z_k')/F_var;
            
            D = (D_1 - (D_2*(D_3\D_4)));
            V_univariate = (v_k*lm.z_k'/F_var)*(pinv(D)*((lm.z_k*v_k')/F_var));
            V_univariate = max(real(V_univariate), 0); % sanitise for numerical rounding errors
        end
    end
end