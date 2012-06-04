% Lagrange Multiplier fMRI model tests example script

% 1) load LM_testData.mat

% 2) reshape data for input
%   data dimensions currently arranged 4D as (x,y,z,time) for x,y,z voxels
%   but this software wants it 2D as (time, voxels) so we reshape as
dShape = size(data);
for t=1:dShape(4)
    dbuf = data(:,:,:,t);
    d(t,:) = dbuf(:).';
end

% 3) Run the model violation tests as e.g.:

% initialise the test object for Double Gamma test
test_DG = LM_DG(d, input_stim_init, TR);  
% adjust any parameters here if you want..

% then run the test for Double Gamma model violation
test_DG = test_DG.apply_test_across_data;   

% and similarly for other tests....
test_NL = LM_NL(d, input_stim_init, TR);  % HRF non-linearity test
test_NL = test_NL.apply_test_across_data;
test_TV = LM_TV(d, input_stim_init, TR);  % time-varying HRF test
test_TV = test_TV.apply_test_across_data;
test_DG_SPM = LM_DG(d, input_stim_init, TR, 'SPM'); 
test_DG_SPM = test_DG_SPM.apply_test_across_data;
% a different way to change test options
test_DG_FSL = test_DG; % copy previous test options
test_DG_FSL.TEST_DERIVATIVE = 'FSL';
test_DG_FSL = test_DG_FSL.apply_test_across_data;   

% now to reconstruct the output test statistics into a spatial map
V_DG = reshape(test_DG.Vt, [80 80 1]);
V_DG(isnan(V_DG)) = -max(V_DG(:));
V_NL = reshape(test_NL.Vt, [80 80 1]);
V_NL(isnan(V_NL)) = -max(V_NL(:));
V_TV = reshape(test_TV.Vt, [80 80 1]);
V_TV(isnan(V_DG)) = -max(V_TV(:));
V_DG_SPM = reshape(test_DG_SPM.Vt, [80 80 1]);
V_DG_SPM(isnan(V_DG_SPM)) = -max(V_DG_SPM(:));
V_DG_FSL = reshape(test_DG_FSL.Vt, [80 80 1]);
V_DG_FSL(isnan(V_DG_FSL)) = -max(V_DG_FSL(:));

% and have a look at the results
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
subplot(2,3,1)
imagesc(data(:,:,1,1));
title('Structural')
axis off
subplot(2,3,2)
imagesc(V_DG)
colorbar
title('Double Gamma (Canonical) HRF anomaly')
axis off
subplot(2,3,3)
imagesc(V_DG_SPM)
colorbar
title('DG + SPM derivative, HRF anomaly')
axis off
subplot(2,3,4)
imagesc(V_DG_FSL)
colorbar
title('DG + FSL derivative, HRF anomaly')
axis off
subplot(2,3,5)
imagesc(V_NL)
colorbar
title('Non-linear HRF anomaly')
axis off
subplot(2,3,6)
imagesc(V_TV)
colorbar
title('Time-varying HRF anomaly')
axis off