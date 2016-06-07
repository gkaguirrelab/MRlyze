function convert_Mathematica_fine_templates(session_dir,tdir)

% Takes the .mgz output from Mathematica, convertes to nii.gz and separates
%   out the pol, ecc, and areas maps into individual volumes.
%
%   Usage: convert_Mathematica_fine_templates(session_dir,tdir)
%
%   Written by Andrew S Bock Jul 2015

%% Find template files
t = listdir(fullfile(tdir,'*.mgz'),'files');
%% Template fit to pRF ROI
disp('Creating pRF template ROIs...');
for dd = 1:length(t)
    hemi = t{dd}(1:2);
    % Load in a temporary file to overwrite
    tmp = load_nifti(fullfile(session_dir,'anat_templates','lh.areas.anat.nii.gz'));
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
end
disp('done.');