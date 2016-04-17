function create_pRF_fine_template(session_dir,subject_name,hemis)
% Takes the .mgz output from Mathematica, convertes to nii.gz and separates
%   out the pol, ecc, and areas maps into individual volumes.
%
%   Usage: create_pRF_template(session_dir,subject_name,hemi)
%
%   Written by Andrew S Bock Apr 2015

%% Set default parameters
if ~exist('session_dir','var')
    error('No ''session_dir'' defined')
end
if ~exist('subject_name','var')
    error('No ''subject_name'' defined')
end
if ~exist('hemi','var')
    hemis = {'lh' 'rh'};
end
%% Add to log
SaveLogInfo(session_dir,mfilename,session_dir,subject_name,hemis)

%% Template fit to pRF ROI
disp('Creating pRF template ROIs...');
for hh = 1:length(hemis)
    % Load in a temporary file to overwrite
    tmp = load_nifti(fullfile(session_dir,'lh.areas.nii.gz'));
    % Load the output from Mathematica
    tmpmgh = load_mgh(fullfile(session_dir,'pRFs',[hemis{hh} '_fine_template_pRFs.mgz']));
    % Pol
    tmp.vol = squeeze(tmpmgh(1,:,1))';
    tmp.vol(tmp.vol==99999) = nan;
    tmp.vol = deg2rad(tmp.vol) - pi/2;
    save_nifti(tmp,fullfile(session_dir,[hemis{hh} '.pol_fine_pRF.nii.gz']));
    % Ecc
    tmp.vol = squeeze(tmpmgh(1,:,2))';
    tmp.vol(tmp.vol==99999) = nan;
    save_nifti(tmp,fullfile(session_dir,[hemis{hh} '.ecc_fine_pRF.nii.gz']));
    % Areas
    tmp.vol = squeeze(tmpmgh(1,:,3))';
    tmp.vol(tmp.vol==99999) = nan;
    save_nifti(tmp,fullfile(session_dir,[hemis{hh} '.areas_fine_pRF.nii.gz']));
    if strcmp(hemis{hh},'rh')
        [~,~] = system(['mri_surf2surf --srcsubject ' subject_name ' --sval ' ...
            fullfile(session_dir,'rh.pol_fine_pRF.nii.gz') ' ' ...
            '--trgsubject ' subject_name '/xhemi --tval ' ...
            fullfile(session_dir,'rh.pol_fine_pRF.nii.gz') ' --hemi lh']);
        [~,~] = system(['mri_surf2surf --srcsubject ' subject_name ' --sval ' ...
            fullfile(session_dir,'rh.ecc_fine_pRF.nii.gz') ' ' ...
            '--trgsubject ' subject_name '/xhemi --tval ' ...
            fullfile(session_dir,'rh.ecc_fine_pRF.nii.gz') ' --hemi lh']);
        [~,~] = system(['mri_surf2surf --srcsubject ' subject_name ' --sval ' ...
            fullfile(session_dir,'rh.areas_fine_pRF.nii.gz') ' ' ...
            '--trgsubject ' subject_name '/xhemi --tval ' ...
            fullfile(session_dir,'rh.areas_fine_pRF.nii.gz') ' --hemi lh']);
        % convert to rh polar angle
        polname = fullfile(session_dir,'rh.pol_fine_pRF.nii.gz');
        tmp = load_nifti(polname);
        upper = tmp.vol>=0;
        lower = tmp.vol<0;
        tmp.vol(upper) = -(tmp.vol(upper) - pi);
        tmp.vol(lower) = -(tmp.vol(lower) + pi);
        save_nifti(tmp,polname);
        % clear up areas
        areaname = fullfile(session_dir,'rh.areas_fine_pRF.nii.gz');
        tmp = load_nifti(areaname);
        tmp.vol = -tmp.vol;
        save_nifti(tmp,areaname);
    end
end
disp('done.');