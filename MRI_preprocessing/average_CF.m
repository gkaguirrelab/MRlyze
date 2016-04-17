function average_CF(session_dir,subject_name,runNums,hemis,srcROI,template,srcfunc,trgfunc,cond,V1only)

%   Average CF maps
%
%   Usage:
%   average_CF(session_dir,subject_name,runNums,hemis,srcROI,template,srcfunc,trgfunc,cond,V1only)
%
%   Written by Andrew S Bock Jan 2016

%% Set defaults
if ~exist('srcROI','var')
    srcROI = 'cortex';
end
if ~exist('template','var')
    template = 'fine';
end
if ~exist('hemis','var')
    hemis = {'lh' 'rh'};
end
if ~exist('V1only','var')
    V1only = 1; % default, use V1. 0 = V1-V3
end
if ~exist('maps','var')
    maps = {'R2' 'shiftt' 'trgpeakt' 'sig1' 'sig2' 'sig3' 'sig4' 'ecc' 'pol'};
end
disp(['Session_dir = ' session_dir]);
disp(['Working on runs ' num2str(runNums)]);
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
%% Average maps
progBar = ProgressBar(length(hemis),'Averaging maps...');
d = find_bold(session_dir);
for hh = 1:length(hemis)
    clear sum_mat
    hemi = hemis{hh};
    for m = 1:length(maps)
        tmp = [];
        % Get the map values for all the runs
        for rr = runNums
            if V1only
                tmpfile = fullfile(session_dir,'CFs',d{rr},...
                    [hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc ...
                    '.' cond '.run' num2str(rr) '.' maps{m} '.V1.cfs']);
            else
                tmpfile = fullfile(session_dir,'CFs',d{rr},...
                    [hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc ...
                    '.' cond '.run' num2str(rr) '.' maps{m} '.cfs']);
            end
            if strcmp(srcROI,'cortex')
                tmpvol = load_nifti([tmpfile '.nii.gz']);
            else
                tmpvol = load_nifti([tmpfile '.orig.nii.gz']);
            end
            tmpvect = reshape(tmpvol.vol,...
                size(tmpvol.vol,1)*size(tmpvol.vol,2)*size(tmpvol.vol,3),1);
            tmp = [tmp,tmpvect];
        end
        sum_mat(m,:,:) = tmp; % create summary matrix of all maps and runs
    end
    % Average the maps across the runs
    for m = 1:length(maps)
        tmp = squeeze(sum_mat(m,:,:));
        tmpavg = [];
        % load template nii
        if strcmp(srcROI,'cortex')
            nii = load_nifti([tmpfile '.nii.gz']);
        else
            nii = load_nifti([tmpfile '.orig.nii.gz']);
        end
        if strcmp('copol',maps{m}) || strcmp('varpol',maps{m})
            %tmpavg = nancirc_mean(tmp);
            tmpavg = circ_mean(tmp');
            tmpavg = tmpavg';
            tmpavg = reshape(tmpavg,size(nii.vol));
            nii.vol = tmpavg;
        else
            %tmpavg = nanmean(tmp,2);
            if strcmp('co',maps{m})
                tmp = fisher_z_corr(tmp);
                tmpavg = mean(tmp,2);
                tmpavg = reshape(tmpavg,size(nii.vol));
                nii.vol = tmpavg;
            else
                tmpavg = mean(tmp,2);
                tmpavg = reshape(tmpavg,size(nii.vol));
                nii.vol = tmpavg;
            end
        end
        % Save nifti
        if V1only
            save_nifti(nii,fullfile(session_dir,'CFs',...
                [hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                cond '.avg.' maps{m} '.V1.cfs.nii.gz']));
        else
            save_nifti(nii,fullfile(session_dir,'CFs',...
                [hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                cond '.avg.' maps{m} '.cfs.nii.gz']));
        end
        % Project to fsaverage_sym space
        if strcmp(srcROI,'cortex')
            if V1only
                sval = fullfile(session_dir,'CFs',...
                    [hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                    cond '.avg.' maps{m} '.V1.cfs.nii.gz']);
                oval = fullfile(session_dir,'CFs',...
                    [hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                    cond '.avg.' maps{m} '.V1.cfs.sym.nii.gz']);
            else
                sval = fullfile(session_dir,'CFs',...
                    [hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                    cond '.avg.' maps{m} '.cfs.nii.gz']);
                oval = fullfile(session_dir,'CFs',...
                    [hemi '.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                    cond '.avg.' maps{m} '.cfs.sym.nii.gz']);
            end
            if strcmp(hemi,'lh')
                system(['mri_surf2surf --hemi lh --srcsubject ' subject_name ...
                    ' --sval ' sval ' --trgsubject fsaverage_sym --trgsurfval ' ...
                    oval]);
            elseif strcmp(hemi,'rh') % need to complete lh first
                system(['mri_surf2surf --hemi lh --srcsubject ' subject_name ...
                    '/xhemi --sval ' sval ' --trgsubject fsaverage_sym --trgsurfval ' ...
                    oval]);
                % average lh and rh
                if V1only
                    lhvol = fullfile(session_dir,'CFs',...
                        ['lh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.V1.cfs.sym.nii.gz']);
                    rhvol = fullfile(session_dir,'CFs',...
                        ['rh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.V1.cfs.sym.nii.gz']);
                else
                    lhvol = fullfile(session_dir,'CFs',...
                        ['lh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.cfs.sym.nii.gz']);
                    rhvol = fullfile(session_dir,'CFs',...
                        ['rh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.cfs.sym.nii.gz']);
                end
                lh = load_nifti(lhvol);
                rh = load_nifti(rhvol);
                if strcmp('copol',maps{m}) || strcmp('varpol',maps{m}) || ...
                        strcmp('pol',maps{m})
                    upper = rh.vol>=0;
                    lower = rh.vol<0;
                    rh.vol(upper) = -(rh.vol(upper) - pi);
                    rh.vol(lower) = -(rh.vol(lower) + pi);
                    tmp = [lh.vol,rh.vol];
                    %tmpavg = nancirc_mean(tmp);
                    tmpavg = circ_mean(tmp');
                    tmpavg = tmpavg';
                    lh.vol = tmpavg;
                else
                    tmp = [lh.vol,rh.vol];
                    %tmpavg = nanmean(tmp,2);
                    tmpavg = mean(tmp,2);
                    lh.vol = tmpavg;
                end
                % Save nifti
                if V1only
                    save_nifti(lh,fullfile(session_dir,'CFs',...
                        ['mh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.V1.cfs.sym.nii.gz']));
                    % Project back to subject space
                    msval = fullfile(session_dir,'CFs',...
                        ['mh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.V1.cfs.sym.nii.gz']);
                    moval = fullfile(session_dir,'CFs',...
                        ['mh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.V1.cfs.nii.gz']);
                    [~,~] = system(['mri_surf2surf --hemi lh --srcsubject fsaverage_sym ' ...
                        ' --sval ' msval ' --trgsubject ' subject_name ' --trgsurfval ' ...
                        moval]);
                else
                    save_nifti(lh,fullfile(session_dir,'CFs',...
                        ['mh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.cfs.sym.nii.gz']));
                    % Project back to subject space
                    msval = fullfile(session_dir,'CFs',...
                        ['mh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.cfs.sym.nii.gz']);
                    moval = fullfile(session_dir,'CFs',...
                        ['mh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.cfs.nii.gz']);
                    [~,~] = system(['mri_surf2surf --hemi lh --srcsubject fsaverage_sym ' ...
                        ' --sval ' msval ' --trgsubject ' subject_name ' --trgsurfval ' ...
                        moval]);
                end
            end
        else
            if strcmp(hemi,'rh') % need to complete lh first
                if V1only
                    lh_vol_in = fullfile(session_dir,'CFs',...
                        ['lh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.V1.cfs.nii.gz']);
                    rh_vol_in = fullfile(session_dir,'CFs',...
                        ['rh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.V1.cfs.nii.gz']);
                    mh_vol_out = fullfile(session_dir,'CFs',...
                        ['mh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.V1.cfs.mni.nii.gz']);
                    average_in_cvsMNI(session_dir,subject_name,maps{m},lh_vol_in,rh_vol_in,mh_vol_out)
                    % Project back to subject space
                    out_vol = fullfile(session_dir,'CFs',...
                        ['mh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.V1.cfs.nii.gz']);
                    ref_vol = fullfile(session_dir,'CFs',...
                        ['lh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.V1.cfs.nii.gz']);
                    apply_cvs_inverse(mh_vol_out,out_vol,ref_vol,subject_name,SUBJECTS_DIR)
                else
                    lh_vol_in = fullfile(session_dir,'CFs',...
                        ['lh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.cfs.nii.gz']);
                    rh_vol_in = fullfile(session_dir,'CFs',...
                        ['rh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.cfs.nii.gz']);
                    mh_vol_out = fullfile(session_dir,'CFs',...
                        ['mh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.cfs.mni.nii.gz']);
                    average_in_cvsMNI(session_dir,subject_name,maps{m},lh_vol_in,rh_vol_in,mh_vol_out)
                    % Project back to subject space
                    out_vol = fullfile(session_dir,'CFs',...
                        ['mh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.cfs.nii.gz']);
                    ref_vol = fullfile(session_dir,'CFs',...
                        ['lh.' srcROI '.' template '.' srcfunc '.' trgfunc '.' ...
                        cond '.avg.' maps{m} '.cfs.nii.gz']);
                    apply_cvs_inverse(mh_vol_out,out_vol,ref_vol,subject_name,SUBJECTS_DIR)
                end
            end
        end
    end
    progBar(hh);
end
disp('done.');