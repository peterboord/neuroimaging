function pesticaStage5(epi,epi_mask,reference,matlabCallNr)
% Usage:
% pesticaStage5(epi,epi_mask,reference,matlabCallNr)
% matlabCallNr = 'testing','1', or '2'
% where 1 & 2 refer to calling instance

% set paths
MATLAB_AFNI_DIR = getenv('MATLAB_AFNI_DIR');
addpath(MATLAB_AFNI_DIR);
PESTICA_DIR = getenv('PESTICA_DIR');
addpath(PESTICA_DIR);
addpath(fullfile(PESTICA_DIR,'matlab_retroicor'));
% use the data to adjust the dithering of each cardiac peak (when using parallel monitored pulse, this shouldn't do much,
% but this does help the PESTICA cardiac estimators)
% adjusting and getting IRFs - second round
load impulse_responses.mat;
disp('Wait, script starting...');
PESTICA_SLICE_TIMING = getenv('PESTICA_SLICE_TIMING');
if strcmp(matlabCallNr,'testing')
    disp(size(cardph));
    disp(size(coeffs_card));
    disp([epi,'+orig']);
    disp([epi_mask,'+orig']);
    disp(size(cmask));
    disp(PESTICA_SLICE_TIMING);
    disp(['../',reference]);
else
    switch matlabCallNr
        case '1'
            fit_each_physio_peak(cardph,coeffs_card,[epi,'+orig'],[epi_mask,'+orig'],cmask,PESTICA_SLICE_TIMING,['../',reference]);
            load physio_fitted.mat;
            retroicor_get_irf([epi,'+orig'],cardphase,respph,5,PESTICA_SLICE_TIMING,[epi_mask,'+orig']);
            disp('Stage 5 halfway done!');
        case '2'
            % Applying IRFs and correcting data with IRF-RETROICOR
            % writing out corrected data as: irfretroicor
            % writing out statistical coupling maps as: coupling_irfret_card coupling_irfret_resp
            irf_retroicor([epi,'+orig'],cardph,respph,coeffs_card,coeffs_resp,PESTICA_SLICE_TIMING,[epi_mask,'+orig'],['../',reference]);
            disp('Stage 5 Done!');
        otherwise
            error('wrong matlabCallNr');
    end
end
