function decimate_templates(subject_name,tdir,src_surf,trg_surf)

% Decimates pRF templates.  Assumes that the following has been done in
% terminal (must be Linux!).
%
%   cd $SUBJECTS_DIR/<subject_name>/surf
%   mris_decimate -d 0.1 ./lh.inflated ./lh.0.1.inflated
%   mris_decimate -d 0.1 ./rh.inflated ./rh.0.1.inflated
%
%   Usage:
%   decimate_templates(subject_name,tdir,src_surf,trg_surf)
%
%   Written by Andrew S Bock Sep 2015

%% set defaults
if ~exist(fullfile(tdir,'decimated_templates'),'dir')
    mkdir(fullfile(tdir,'decimated_templates'));
end
templates = listdir(fullfile(tdir,'*.nii.gz'),'files');
if ~exist('src_surf','var')
    src_surf = 'inflated';
end
if ~exist('trg_surf','var')
    trg_surf = '0.1.inflated';
end
lhct = 0;
rhct = 0;
%%
disp('Decimating templates...');
for i = 1:length(templates)
    hemi = templates{i}(1:2);
    if strcmp(hemi,'lh');
        lhct = lhct + 1;
    elseif strcmp(hemi,'rh')
        rhct = rhct + 1;
    end
    if lhct == 1 || rhct == 1
        [sind] = find_closest_verts(subject_name,hemi,src_surf,trg_surf);
    end
    in_name = fullfile(tdir,templates{i});
    out_name = fullfile(tdir,'decimated_templates',templates{i});
    project_2_decimate(in_name,out_name,sind);
end
disp('done.');