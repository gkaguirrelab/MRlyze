function do_pRF(session_dir,subject_name,runNum,hemi,srcROI,func,srcind,srcfile,fieldSize,imFile,paramsFile)

% Takes in denoised timeseries, fits pRF using retinotopic mapping stimuli.
%   The main function is 'pRF.m'
%
%   Usage: do_pRF(session_dir,subject_name,runs,srcROI)
%
% The 'runs' input is a vector specifying the runs (i.e. bold directories
%   in the session_dir) which contain retinotopic mapping data.
%
%   e.g. runs = [1 3 5];
%
% Defaults:
%   srcROI - can either be 'occipital', 'cortex' or 'subcortical'
%
% Assumes data files containing the images and parameters used for
% retinotopic mapping are stored as follows:
%
%   e.g. for run3
%       images - fullfile(session_dir,'Stimuli','run3','bars_images.mat');
%       params - fullfile(session_dir,'Stimuli','run3','bars_params.mat');
%
%
% Output maps can be found in the 'out_dir', which is set as:
%
%   out_dir = fullfile(session_dir,'prfs')
%
%   Written by Andrew S Bock Jan 2015

%% Set defaults
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
%% Add to log
SaveLogInfo(session_dir,mfilename,session_dir,subject_name,runNum,hemi,srcROI,srcind,srcfile,fieldSize,imFile,paramsFile);

%% Find bold run directories
d = find_bold(session_dir);
disp(['Session_dir = ' session_dir]);
%% Load timecourses
disp('Loading timecourses...');
% load files
src = load_nifti(srcfile);
srcdims = size(src.vol);
nFrames = srcdims(4);
TR = src.pixdim(5)/1000;
if TR < 0.1
   error('TR is less than 0.1, header TR most likely not in msec'); 
end
srctc = reshape(src.vol,srcdims(1)*srcdims(2)*srcdims(3),srcdims(4))';
% Pull out relevant timecourses
srctc = srctc(:,srcind);
% Set timecourses with very little variation (var<.1) to flat
srctc = set_to_flat(srctc);
disp('done.');
%% Run pRF
[prfs] = calc_pRF(srctc,TR,nFrames,fieldSize,imFile,paramsFile);
% Save pRF data
if ~exist(fullfile(session_dir,'pRFs'),'dir')
    mkdir(fullfile(session_dir,'pRFs'));
end
if ~exist(fullfile(session_dir,'pRFs',d{runNum}),'dir')
    mkdir(fullfile(session_dir,'pRFs',d{runNum}));
end
save(fullfile(session_dir,'pRFs',d{runNum},[hemi '.' srcROI '.' func '.prfs.mat']),'prfs');
%% Plot in src volume
disp('Making pRF maps...');
if strcmp(srcROI,'cortex');
    mri = load_nifti(fullfile(session_dir,[hemi '.ecc.nii.gz']));
else
    mri = load_nifti(fullfile(session_dir,d{runNum},'single_TR.nii.gz'));
end
% Save co and sig maps
mri.vol = nan(size(mri.vol));
co = mri;
copeakt = mri;
varr = mri;
varpeakt = mri;
var_as_co = mri;
cosig1 = mri;
cosig2 = mri;
cosig3 = mri;
cosig4 = mri;
varsig1 = mri;
varsig2 = mri;
varsig3 = mri;
varsig4 = mri;
disp('Saving correlation and variance explained maps...');
%co and copeakt
co.vol(srcind) = prfs.co;
save_nifti(co,fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.co.prfs.nii.gz']));
copeakt.vol(srcind) = prfs.copeakt;
save_nifti(copeakt,fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.copeakt.prfs.nii.gz']));
%var and varpeakt
varr.vol(srcind) = prfs.var;
save_nifti(varr,fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.var.prfs.nii.gz']));
varpeakt.vol(srcind) = prfs.varpeakt;
save_nifti(varpeakt,fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.varpeakt.prfs.nii.gz']));
%sig1
cosig1.vol(srcind) = prfs.cosig(:,1);
save_nifti(cosig1,fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.cosig1.prfs.nii.gz']));
%sig2
cosig2.vol(srcind) = prfs.cosig(:,2);
save_nifti(cosig2,fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.cosig2.prfs.nii.gz']));
%sig3
cosig3.vol(srcind) = prfs.cosig(:,3);
save_nifti(cosig3,fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.cosig3.prfs.nii.gz']));
%sig4
cosig4.vol(srcind) = prfs.cosig(:,4);
save_nifti(cosig4,fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.cosig4.prfs.nii.gz']));
% For variance, flip Gaussian params so all correlations are positive
% Find var_as_co values that are negative (these we will flip)
flip_sigs_ind = prfs.var_as_co;
flip_sigs_ind(isnan(flip_sigs_ind)) = 0;
flip_sigs_ind(flip_sigs_ind>0) = 0;
flip_sigs_ind(flip_sigs_ind<0) = 1;
flip_sigs_ind = logical(flip_sigs_ind);
% Create the flipped values (i.e. if positive correlation)
newsig3 = -prfs.varsig(flip_sigs_ind,3);
newsig4 = -prfs.varsig(flip_sigs_ind,4);
newvar_as_co = -prfs.var_as_co(flip_sigs_ind);
% Replace with flipped values
prfs.varsig(flip_sigs_ind,3) = newsig3;
prfs.varsig(flip_sigs_ind,4) = newsig4;
prfs.var_as_co(flip_sigs_ind) = newvar_as_co;
% Save Gaussian params
disp('Saving Gaussian parameter maps...');
var_as_co.vol(srcind) = prfs.var_as_co;
save_nifti(var_as_co,fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.var_as_co.prfs.nii.gz']));
%sig1
varsig1.vol(srcind) = prfs.varsig(:,1);
save_nifti(varsig1,fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.varsig1.prfs.nii.gz']));
%sig2
varsig2.vol(srcind) = prfs.varsig(:,2);
save_nifti(varsig2,fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.varsig2.prfs.nii.gz']));
%sig3
varsig3.vol(srcind) = prfs.varsig(:,3);
save_nifti(varsig3,fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.varsig3.prfs.nii.gz']));
%sig4
varsig4.vol(srcind) = prfs.varsig(:,4);
save_nifti(varsig4,fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.varsig4.prfs.nii.gz']));
% Save pol and ecc maps
disp('Saving eccentricity and polar angle maps...');
mri.vol = nan(size(mri.vol));
coecc = mri;
copol = mri;
varecc = mri;
varpol = mri;
coecc.vol(srcind) = prfs.coecc;
save_nifti(coecc,fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.coecc.prfs.nii.gz']));
copol.vol(srcind) = prfs.copol;
save_nifti(copol,fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.copol.prfs.nii.gz']));
varecc.vol(srcind) = prfs.varecc;
save_nifti(varecc,fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.varecc.prfs.nii.gz']));
varpol.vol(srcind) = prfs.varpol;
save_nifti(varpol,fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.varpol.prfs.nii.gz']));
if ~strcmp(srcROI,'cortex');
    % Project to anatomical space
    disp('Projecting maps to anatomical space...');
    maps = {'co' 'copeakt' 'cosig1' 'cosig2' 'cosig3' 'cosig4' 'coecc' 'copol' ...
            'var' 'varpeakt' 'var_as_co' 'varsig1' 'varsig2' 'varsig3' 'varsig4' 'varecc' 'varpol'};
    for m = 1:length(maps);
        savename = fullfile(session_dir,[hemi '.' srcROI '.run' num2str(runNum) '.' maps{m} '.prfs']);
        [~,~] = system(['mri_vol2vol --mov ' savename '.nii.gz ' ...
            '--targ ' fullfile(SUBJECTS_DIR,subject_name,'mri','orig.mgz') ...
            ' --o ' savename '.orig.nii.gz --reg ' ...
            fullfile(session_dir,d{runNum},'brf_bbreg.dat') ...
            ' --interp nearest']);
    end
end
disp('done.');