function convert_Mathematica_template(session_dir,template_dir,tNum)

% Takes the .mgz output from Mathematica, convertes to nii.gz and separates
%   out the pol, ecc, and areas maps into individual volumes.
%
%   Usage: convert_Mathematica_templates(session_dir,template_dir)
%
%   Written by Andrew S Bock Jul 2015

%% Set default parameters
if ~exist('template_dir','var')
    template_dir = fullfile(session_dir,'pRFs','coarse_model_templates');
end
%% Find template files
t = listdir(fullfile(template_dir,'*.mgz'),'files');
%% Template fit to pRF ROI
disp('Creating templates...');
for dd = tNum
    hemi = t{dd}(1:2);
    % Load in a temporary file to overwrite
    inTemp = fullfile(session_dir,'anat_templates',[hemi '.areas.anat.nii.gz']);
    outDir  = fullfile(session_dir,'anat_templates',sprintf('temp%05d',dd));
    mkdir(outDir);
    outTemp = fullfile(outDir,sprintf('temp%05d.nii.gz',dd));
    system(['cp ' inTemp ' ' outTemp]);
    pause(5);
    tmp = load_nifti(outTemp);
    % Load the output from Mathematica
    tmpmgh = load_mgh(fullfile(template_dir,t{dd}));
    % get template name
    endname = strfind(t{dd},'.mgz');
    fname = t{dd}(4:endname-1);
    % Areas
    Areas = tmp;
    Areas.vol = squeeze(tmpmgh(1,:,3))';
    Areas.vol(Areas.vol==99999) = nan;
    Areas.vol(abs(Areas.vol)>3) = nan;
    save_nifti(Areas,fullfile(template_dir,[hemi '.areas.' fname '.nii.gz']));
    % Pol
    Pol = tmp;
    Pol.vol = squeeze(tmpmgh(1,:,1))';
    Pol.vol(Pol.vol==99999) = nan;
    Pol.vol(abs(Areas.vol)>3) = nan;
    Pol.vol = deg2rad(Pol.vol) - pi/2;
    save_nifti(Pol,fullfile(template_dir,[hemi '.pol.' fname '.nii.gz']));
    % Ecc
    Ecc = tmp;
    Ecc.vol = squeeze(tmpmgh(1,:,2))';
    Ecc.vol(Ecc.vol==99999) = nan;
    Ecc.vol(abs(Areas.vol)>3) = nan;
    save_nifti(Ecc,fullfile(template_dir,[hemi '.ecc.' fname '.nii.gz']));
    % RIGHT HEMISPHERE
    if strcmp(hemi,'rh')
        % convert to rh polar angle
        polname = fullfile(template_dir,[hemi '.pol.' fname '.nii.gz']);
        tmp = load_nifti(polname);
        upper = tmp.vol>=0;
        lower = tmp.vol<0;
        tmp.vol(upper) = (tmp.vol(upper) - pi);
        tmp.vol(lower) = (tmp.vol(lower) + pi);
        save_nifti(tmp,polname);
        % fix rh areas
        areaname = fullfile(template_dir,[hemi '.areas.' fname '.nii.gz']);
        tmp = load_nifti(areaname);
        tmp.vol = -tmp.vol;
        save_nifti(tmp,areaname);
    end
    system(['rm -rf ' outDir]);
end
disp('done.');