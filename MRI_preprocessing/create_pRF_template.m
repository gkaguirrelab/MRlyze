function create_pRF_template(session_dir,subject_name,hemis)
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
outDir = fullfile(session_dir,'pRFs','pRF_templates');
if ~exist(outDir,'dir')
    mkdir(outDir);
end
for hh = 1:length(hemis)
    hemi = hemis{hh};
    % Load in a temporary file to overwrite
    tmp = load_nifti(fullfile(session_dir,'pRFs','anat_templates','lh.areas.anat.nii.gz'));
    % Load the output from Mathematica
    tmpmgh = load_mgh(fullfile(session_dir,'pRFs',[hemi '_template_pRFs.mgz']));
    % Pol
    tmp.vol = squeeze(tmpmgh(1,:,1))';
    tmp.vol(tmp.vol==99999) = nan;
    tmp.vol = deg2rad(tmp.vol) - pi/2;
    save_nifti(tmp,fullfile(outDir,[hemi '.pol.pRF.nii.gz']));
    % Ecc
    tmp.vol = squeeze(tmpmgh(1,:,2))';
    tmp.vol(tmp.vol==99999) = nan;
    save_nifti(tmp,fullfile(outDir,[hemi '.ecc.pRF.nii.gz']));
    % Areas
    tmp.vol = squeeze(tmpmgh(1,:,3))';
    tmp.vol(tmp.vol==99999) = nan;
    save_nifti(tmp,fullfile(outDir,[hemi '.areas.pRF.nii.gz']));
    if strcmp(hemi,'rh')
        %         [~,~] = system(['mri_surf2surf --srcsubject ' subject_name ' --sval ' ...
        %             fullfile(session_dir,'rh.pol_pRF.nii.gz') ' ' ...
        %             '--trgsubject ' subject_name '/xhemi --tval ' ...
        %             fullfile(session_dir,'rh.pol_pRF.nii.gz') ' --hemi lh']);
        %         [~,~] = system(['mri_surf2surf --srcsubject ' subject_name ' --sval ' ...
        %             fullfile(session_dir,'rh.ecc_pRF.nii.gz') ' ' ...
        %             '--trgsubject ' subject_name '/xhemi --tval ' ...
        %             fullfile(session_dir,'rh.ecc_pRF.nii.gz') ' --hemi lh']);
        %         [~,~] = system(['mri_surf2surf --srcsubject ' subject_name ' --sval ' ...
        %             fullfile(session_dir,'rh.areas_pRF.nii.gz') ' ' ...
        %             '--trgsubject ' subject_name '/xhemi --tval ' ...
        %             fullfile(session_dir,'rh.areas_pRF.nii.gz') ' --hemi lh']);
        % convert to rh polar angle
        polname = fullfile(outDir,[hemi '.pol.pRF.nii.gz']);
        tmp = load_nifti(polname);
        upper = tmp.vol>=0;
        lower = tmp.vol<0;
        tmp.vol(upper) = (tmp.vol(upper) - pi);
        tmp.vol(lower) = (tmp.vol(lower) + pi);
        save_nifti(tmp,polname);
        % fix rh areas
        areaname = fullfile(outDir,[hemi '.areas.pRF.nii.gz']);
        tmp = load_nifti(areaname);
        tmp.vol = -tmp.vol;
        save_nifti(tmp,areaname);
    end
end
disp('done.');