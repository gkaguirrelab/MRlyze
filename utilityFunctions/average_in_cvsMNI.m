function average_in_cvsMNI(session_dir,subject_name,map_type,lh_vol_in,rh_vol_in,mh_vol_out)

% Note: requires a copy of subject called "<subject_name>_xhemi".  The way
% that the cvs registration works is not compatible with the normal /xhemi
% directory.
%
%   Written by Andrew S Bock May 2015

%% Set temporary output files
[~,lhname,lhext] = fileparts(lh_vol_in);
if strcmp(lhext,'.gz')
    [~,lhname,~] = fileparts(lhname);
    lh_vol_out = fullfile(session_dir,[lhname '_tmp.nii.gz']);
elseif strcmp(lhext,'.nii')
    lh_vol_out = fullfile(session_dir,[lhname '_tmp.nii.gz']);
end
[~,rhname,rhext] = fileparts(rh_vol_in);
if strcmp(rhext,'.gz')
    [~,rhname,~] = fileparts(rhname);
    rh_vol_out = fullfile(session_dir,[rhname '_tmp.nii.gz']);
elseif strcmp(rhext,'.nii')
    rh_vol_out = fullfile(session_dir,[rhname '_tmp.nii.gz']);
end
SUBJECTS_DIR=getenv('SUBJECTS_DIR');
%% Average hemispheres in MNI space
% left hemisphere
apply_cvs(lh_vol_in,lh_vol_out,subject_name);
% Right hemisphere
system(['fslswapdim ' rh_vol_in ' -x y z ' rh_vol_out]);
apply_cvs(rh_vol_out,rh_vol_out,[subject_name '_xhemi']);
% Average hemispheres
disp(['averaging maps into ' mh_vol_out '...']);
lhtmpvol = load_nifti(lh_vol_out);
rhtmpvol = load_nifti(rh_vol_out);
lhtmpvect = reshape(lhtmpvol.vol,size(lhtmpvol.vol,1)*size(lhtmpvol.vol,2)*size(lhtmpvol.vol,3),1);
rhtmpvect = reshape(rhtmpvol.vol,size(rhtmpvol.vol,1)*size(rhtmpvol.vol,2)*size(rhtmpvol.vol,3),1);
tmp = [lhtmpvect,rhtmpvect];
nii = load_nifti(lh_vol_out);
if strcmp('copol',map_type) || strcmp('varpol',map_type)
    %tmpavg = nancirc_mean(tmp);
    tmpavg = circ_mean(tmp');
    tmpavg = tmpavg';
    tmpavg = reshape(tmpavg,size(nii.vol));
    nii.vol = tmpavg;
else
    %tmpavg = nanmean(tmp,2);
    tmpavg = mean(tmp,2);
    tmpavg = reshape(tmpavg,size(nii.vol));
    nii.vol = tmpavg;
end
save_nifti(nii,mh_vol_out);
% Remove temporary output files
delete([lh_vol_out '*']);
delete([rh_vol_out '*']);
disp('done.');