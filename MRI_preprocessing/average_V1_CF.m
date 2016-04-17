function average_V1_CF(session_dir,subject_name,map_type,runs,srcROI,template,hemis)

%   Average CF maps
%
%   Usage:
%   average_CF(session_dir,subject_name,map_type,runs,srcROI,template,hemis)
%
%   e.g.
%   map_type = 'bar'; % 'movie'
%   runs = [1,3,5];
%   average_CF(session_dir,subject_name,map_type,runs)
%
%   defaults:
%   runs = all runs (determined by the number of bold directories.
%   Alternatively, input is a variable runs, e.g. runs = [1,3,5];
%   srcROI = 'cortex';
%   template = 'fine';
%   hemi = {'lh' 'rh'};
%
%   Written by Andrew S Bock Oct 2014

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
if ~exist('maps','var')
    maps = {'R2' 'shiftt' 'trgpeakt' 'sig1' 'sig2' 'sig3' 'sig4' 'ecc' 'pol'};
end
disp(['Session_dir = ' session_dir]);
disp(['Working on runs ' num2str(runs)]);
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
%% Average maps
progBar = ProgressBar(length(hemis),'Averaging maps...');
for hh = 1:length(hemis)
    clear sum_mat
    hemi = hemis{hh};
    for m = 1:length(maps)
        tmp = [];
        % Get the map values for all the runs
        for rr = runs
            tmpfile = fullfile(session_dir,...
                [hemi '.' srcROI '.' template '.run' num2str(rr) '.' maps{m} '.V1cfs']);
            switch srcROI
                case 'cortex'
                    tmpvol = load_nifti([tmpfile '.nii.gz']);
                case 'volume'
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
        switch srcROI
            case 'cortex'
                nii = load_nifti([tmpfile '.nii.gz']);
            case 'volume'
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
        save_nifti(nii,fullfile(session_dir,...
            [hemi '.' srcROI '.' template '.' map_type '.avg.' maps{m} '.V1cfs.nii.gz']));
        % Project to fsaverage_sym space
        switch srcROI
            case 'cortex'
                sval = fullfile(session_dir,...
                    [hemi '.' srcROI '.' template '.' map_type '.avg.' maps{m} '.V1cfs.nii.gz']);
                oval = fullfile(session_dir,...
                    [hemi '.' srcROI '.' template '.' map_type '.avg.' maps{m} '.V1cfs.sym.nii.gz']);
                if strcmp(hemi,'lh')
                    [~,~] = system(['mri_surf2surf --hemi lh --srcsubject ' subject_name ...
                        ' --sval ' sval ' --trgsubject fsaverage_sym --trgsurfval ' ...
                        oval]);
                elseif strcmp(hemi,'rh') % need to complete lh first
                    [~,~] = system(['mri_surf2surf --hemi lh --srcsubject ' subject_name ...
                        '/xhemi --sval ' sval ' --trgsubject fsaverage_sym --trgsurfval ' ...
                        oval]);
                    % average lh and rh
                    lhvol = fullfile(session_dir,...
                        ['lh.' srcROI '.' template '.' map_type '.avg.' maps{m} '.V1cfs.sym.nii.gz']);
                    rhvol = fullfile(session_dir,...
                        ['rh.' srcROI '.' template '.' map_type '.avg.' maps{m} '.V1cfs.sym.nii.gz']);
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
                    save_nifti(lh,fullfile(session_dir,...
                        ['mh.' srcROI '.' template '.' map_type '.avg.' maps{m} '.V1cfs.sym.nii.gz']));
                    % Project back to subject space
                    msval = fullfile(session_dir,...
                        ['mh.' srcROI '.' template '.' map_type '.avg.' maps{m} '.V1cfs.sym.nii.gz']);
                    moval = fullfile(session_dir,...
                        ['mh.' srcROI '.' template '.' map_type '.avg.' maps{m} '.V1cfs.nii.gz']);
                    [~,~] = system(['mri_surf2surf --hemi lh --srcsubject fsaverage_sym ' ...
                        ' --sval ' msval ' --trgsubject ' subject_name ' --trgsurfval ' ...
                        moval]);
                end
            case 'volume'
                if strcmp(hemi,'rh') % need to complete lh first
                    lh_vol_in = fullfile(session_dir,['lh.' srcROI '.' template '.' map_type '.avg.' maps{m} '.V1cfs.nii.gz']);
                    rh_vol_in = fullfile(session_dir,['rh.' srcROI '.' template '.' map_type '.avg.' maps{m} '.V1cfs.nii.gz']);
                    mh_vol_out = fullfile(session_dir,['mh.' srcROI '.' template '.' map_type '.avg.' maps{m} '.V1cfs.mni.nii.gz']);
                    average_in_cvsMNI(session_dir,subject_name,maps{m},lh_vol_in,rh_vol_in,mh_vol_out)
                    % Project back to subject space
                    out_vol = fullfile(session_dir,['mh.' srcROI '.' template '.' map_type '.avg.' maps{m} '.V1cfs.nii.gz']);
                    ref_vol = fullfile(session_dir,['lh.' srcROI '.' template '.' map_type '.avg.' maps{m} '.V1cfs.nii.gz']);
                    system(['mri_vol2vol --inv-morph --targ ' mh_vol_out ...
                        ' --m3z final_CVSmorph_tocvs_avg35_inMNI152.m3z' ...
                        ' --reg ' fullfile(SUBJECTS_DIR,subject_name,'cvs/subj2MNI.dat') ...
                        ' --mov ' ref_vol ' --o ' out_vol]);
                end
        end
    end
    progBar(hh);
end
disp('done.');