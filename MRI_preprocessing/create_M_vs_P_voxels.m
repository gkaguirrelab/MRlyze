function create_M_vs_P_voxels(session_dir,subject_name,func,gfeatNum,hemi,ROI)

%   This will create a binary volume for M and P, for a given ROI
%
%   Usage:
%   create_M_vs_P_voxels(session_dir,subject_name,func,hemi,ROI)
%
%   Written by Andrew S Bock Nov 2015

%% Find bold directories
b = find_bold(session_dir);
featDir = fullfile(session_dir,b{gfeatNum},[func '.gfeat']);
%% Load the M vs P volume
%mp = load_nifti(fullfile(featDir,'cope3.feat','stats','cope1.psc.nii.gz'));
mp = load_nifti(fullfile(featDir,'cope3.feat','stats','zstat1.nii.gz'));

%% Load ROI
switch ROI
    case 'LGN'
        MPfrac = 3; % 1/3 top/bottom
        if strcmp(func(1:2),'s2') || strcmp(func(1:2),'s5')
            area = load_nifti(fullfile(session_dir,[hemi '.s2.' func(4:end) '.LGN.nii.gz']));
        else
            area = load_nifti(fullfile(session_dir,[hemi '.s2.' func '.LGN.nii.gz']));
        end
        MInd = find(area.vol>0);
        PInd = MInd;
    case 'V1'
        MPfrac = 3; % 1/3 of each ROI
        area = load_nifti(fullfile(session_dir,[hemi '.areas.vol.nii.gz']));
        ecc = load_nifti(fullfile(session_dir,[hemi '.ecc.vol.nii.gz']));
        MInd = find((area.vol == 1 | area.vol == -1) & (ecc.vol>10 & ecc.vol<40));
        PInd = find((area.vol == 1 | area.vol == -1) & (ecc.vol<10));
    case 'V2'
        MPfrac = 3; % 1/3 of each ROI
        area = load_nifti(fullfile(session_dir,[hemi '.areas.vol.nii.gz']));
        ecc = load_nifti(fullfile(session_dir,[hemi '.ecc.vol.nii.gz']));
        MInd = find((area.vol == 2 | area.vol == -2) & (ecc.vol>10 & ecc.vol<40));
        PInd = find((area.vol == 2 | area.vol == -2) & (ecc.vol<10));
    case 'V3'
        MPfrac = 3; % 1/3 of each ROI
        area = load_nifti(fullfile(session_dir,[hemi '.areas.vol.nii.gz']));
        ecc = load_nifti(fullfile(session_dir,[hemi '.ecc.vol.nii.gz']));
        MInd = find((area.vol == 3 | area.vol == -3) & (ecc.vol>10 & ecc.vol<40));
        PInd = find((area.vol == 3 | area.vol == -3) & (ecc.vol<10));
    case 'V1_pRF'
        MPfrac = 3; % 1/3 of each ROI
        area = load_nifti(fullfile(session_dir,[hemi '.areas_pRF.vol.nii.gz']));
        ecc = load_nifti(fullfile(session_dir,[hemi '.ecc_pRF.vol.nii.gz']));
        MInd = find((area.vol == 1 | area.vol == -1) & (ecc.vol>10 & ecc.vol<40));
        PInd = find((area.vol == 1 | area.vol == -1) & (ecc.vol<10));
    case 'V2_pRF'
        MPfrac = 3; % 1/3 of each ROI
        area = load_nifti(fullfile(session_dir,[hemi '.areas_pRF.vol.nii.gz']));
        ecc = load_nifti(fullfile(session_dir,[hemi '.ecc_pRF.vol.nii.gz']));
        MInd = find((area.vol == 2 | area.vol == -2) & (ecc.vol>10 & ecc.vol<40));
        PInd = find((area.vol == 2 | area.vol == -2) & (ecc.vol<10));
    case 'V3_pRF'
        MPfrac = 3; % 1/3 of each ROI
        area = load_nifti(fullfile(session_dir,[hemi '.areas_pRF.vol.nii.gz']));
        ecc = load_nifti(fullfile(session_dir,[hemi '.ecc_pRF.vol.nii.gz']));
        MInd = find((area.vol == 3 | area.vol == -3) & (ecc.vol>10 & ecc.vol<40));
        PInd = find((area.vol == 3 | area.vol == -3) & (ecc.vol<10));
end
%% Sort ROI voxels into top (M) and bottom (P) third
% M
Mvals = mp.vol(MInd);
[MY,MI] = sort(Mvals);
fracSize = floor(length(Mvals)/MPfrac);
M.vals = MY(end-(fracSize-1):end);
M.I = MI(end-(fracSize-1):end);
% P
Pvals = mp.vol(PInd);
[PY,PI] = sort(Pvals);
fracSize = floor(length(Pvals)/MPfrac);
P.vals = PY(1:fracSize);
P.I = PI(1:fracSize);
%% Save volumes
tmpM = area;
tmpP = area;
tmpAll = area;
tmpM.vol = zeros(size(tmpM.vol));
tmpP.vol = zeros(size(tmpP.vol));
tmpAll.vol = zeros(size(tmpP.vol));
tmpM.vol(MInd(M.I)) = 1;
tmpP.vol(PInd(P.I)) = 1;
tmpAll.vol([MInd;PInd]) = 1;
save_nifti(tmpM,fullfile(session_dir,[hemi '.' func '.' ROI '.M.nii.gz']));
save_nifti(tmpP,fullfile(session_dir,[hemi '.' func '.' ROI '.P.nii.gz']));
save_nifti(tmpP,fullfile(session_dir,[hemi '.' func '.' ROI '.All.nii.gz']));
tmpM.vol(MInd(M.I)) = mp.vol(MInd(M.I));
tmpP.vol(PInd(P.I)) = mp.vol(PInd(P.I));
tmpAll.vol([MInd;PInd]) = mp.vol([MInd;PInd]);
save_nifti(tmpM,fullfile(session_dir,[hemi '.' func '.' ROI '.M.zstat.nii.gz']));
save_nifti(tmpP,fullfile(session_dir,[hemi '.' func '.' ROI '.P.zstat.nii.gz']));
save_nifti(tmpAll,fullfile(session_dir,[hemi '.' func '.' ROI '.ALL.zstat.nii.gz']));
%% Project to surface, if cortical ROI
inM = fullfile(session_dir,[hemi '.' func '.' ROI '.M.nii.gz']);
inP = fullfile(session_dir,[hemi '.' func '.' ROI '.P.nii.gz']);
inMP = fullfile(session_dir,[func '.M_vs_P.nii.gz']);
outM = fullfile(session_dir,[hemi '.' func '.' ROI '.M.surf.nii.gz']);
outP = fullfile(session_dir,[hemi '.' func '.' ROI '.P.surf.nii.gz']);
outMP = fullfile(session_dir,[hemi '.' func '.' ROI '.M_vs_P.surf.nii.gz']);
if strcmp(ROI,'V1') || strcmp(ROI,'V2') || strcmp(ROI,'V3') || ...
        strcmp(ROI,'V1_pRF') || strcmp(ROI,'V2_pRF') || strcmp(ROI,'V3_pRF')
    system(['mri_vol2surf --src ' inM ...
        ' --regheader ' subject_name ' --hemi ' hemi ...
        ' --out ' outM ' --projfrac 0.5']);
    system(['mri_vol2surf --src ' inP ...
        ' --regheader ' subject_name ' --hemi ' hemi ...
        ' --out ' outP ' --projfrac 0.5']);
    system(['mri_vol2surf --src ' inMP ...
        ' --regheader ' subject_name ' --hemi ' hemi ...
        ' --out ' outMP ' --projfrac 0.5']);
end