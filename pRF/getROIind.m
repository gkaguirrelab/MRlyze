function [ROIind] = getROIind(hemis,ROIs)

for hh = 1:length(hemis)
    hemi = hemis{hh};
    for jj = 1:length(ROIs)
        ROI = ROIs{jj};
        % Get ROIind
 areas = load_nifti(fullfile(session_dir,'anat_templates',[hemi '.areas.anat.vol.nii.gz']));
        ecc = load_nifti(fullfile(session_dir,'anat_templates',[hemi '.ecc.anat.vol.nii.gz']));
        lgn = load_nifti(fullfile(session_dir,'anat_templates',[hemi '.LGN.nii.gz']));
        switch ROI
            case 'V1all'
                ROIind = find(abs(areas.vol)==1); % all of V1
            case 'V1low'
                ROIind = find(abs(areas.vol)==1 & ecc.vol<=5);
            case 'V1mid'
                ROIind = find(abs(areas.vol)==1 & (ecc.vol>5 & ecc.vol<=15));
            case 'V1high'
                ROIind = find(abs(areas.vol)==1 & (ecc.vol>15 & ecc.vol<=40));
            case 'V2andV3'
                ROIind = find(abs(areas.vol)==2 | abs(areas.vol)==3);
            case 'LGN'
                  ROIind = find((lgn.vol)>0);
        end
    end
end