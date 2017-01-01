function savepRF(matFile,outName,outDir,srcROI,srcind,templateFile,subject_name,regFile)

% Calcuate the pRF for a given set of fMRI voxels/vertices
%
%   Usage:
%   savepRF(matFile,outName,outDir,srcROI,srcind,templateFile,subject_name,regFile)
%
%   Note:
%   the variables 'subject_name' and 'regFile' are only relevant for pRFs
%   calculated in the volume
%
%   Written by Andrew S Bock May 2016

%% Set defaults
load(matFile);
SUBJECTS_DIR = getenv('SUBJECTS_DIR');

%% save pRFs in source volume
disp('Saving pRF maps...');
% Load a template file that is in the same space as the desired outputs
if strcmp(srcROI,'cortex');
    mri = load_nifti(templateFile); % e.g. fullfile(session_dir,'anat_templates',[hemi '.ecc.anat.nii.gz'])
else
    mri = load_nifti(templateFile); % e.g. fullfile(session_dir,d{runNum},'single_TR.nii.gz')
end
% Save co and sig maps
mri.vol = nan(size(mri.vol));
co = mri;
copeakt = mri;
cosig1 = mri;
cosig2 = mri;
cosig3 = mri;
cosig4 = mri;
%co and copeakt
co.vol(srcind) = prfs.co;
save_nifti(co,fullfile(outDir,[outName '.' srcROI '.co.prfs.nii.gz']));
copeakt.vol(srcind) = prfs.copeakt;
save_nifti(copeakt,fullfile(outDir,[outName '.' srcROI '.copeakt.prfs.nii.gz']));
%sig1
cosig1.vol(srcind) = prfs.cosig(:,1);
save_nifti(cosig1,fullfile(outDir,[outName '.' srcROI '.cosig1.prfs.nii.gz']));
%sig2
cosig2.vol(srcind) = prfs.cosig(:,2);
save_nifti(cosig2,fullfile(outDir,[outName '.' srcROI '.cosig2.prfs.nii.gz']));
%sig3
cosig3.vol(srcind) = prfs.cosig(:,3);
save_nifti(cosig3,fullfile(outDir,[outName '.' srcROI '.cosig3.prfs.nii.gz']));
%sig4
cosig4.vol(srcind) = prfs.cosig(:,4);
save_nifti(cosig4,fullfile(outDir,[outName '.' srcROI '.cosig4.prfs.nii.gz']));
% Save pol and ecc maps
mri.vol = nan(size(mri.vol));
coecc = mri;
copol = mri;
coecc.vol(srcind) = prfs.coecc;
save_nifti(coecc,fullfile(outDir,[outName '.' srcROI '.coecc.prfs.nii.gz']));
copol.vol(srcind) = prfs.copol;
save_nifti(copol,fullfile(outDir,[outName '.' srcROI '.copol.prfs.nii.gz']));
if ~strcmp(srcROI,'cortex');
    % Project to anatomical space
    disp('Projecting maps to anatomical space...');
    maps = {'co' 'copeakt' 'cosig1' 'cosig2' 'cosig3' 'cosig4' 'coecc' 'copol'};
    for m = 1:length(maps);
        savename = fullfile(outDir,[outName '.' srcROI '.' maps{m} '.prfs']);
        system(['mri_vol2vol --mov ' savename '.nii.gz ' ...
            '--targ ' fullfile(SUBJECTS_DIR,subject_name,'mri','orig.mgz') ...
            ' --o ' savename '.orig.nii.gz --reg ' ...
            regFile ' --interp nearest']);
    end
end
disp('done.');