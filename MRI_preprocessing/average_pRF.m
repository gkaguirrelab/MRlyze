function average_pRF(session_dir,subject_name,runs,srcROI,hemis)

%   Average CF maps
%
%   Usage:
%   average_pRF(session_dir,subject_name,runs,srcROI,hemis)
%
%   e.g.
%   runs = [1,3,5];
%   average_pRF(session_dir,subject_name,runs)
%
%   defaults:
%   runs = all runs (determined by the number of bold directories.
%   Alternatively, input is a variable runs, e.g. runs = [1,3,5];
%   srcROI = 'volume';
%   hemi = {'lh' 'rh'};
%
%   Written by Andrew S Bock Oct 2014

%% Find bold run directories
d = find_bold(session_dir);
nruns = length(d);
if ~exist('runs','var')
    runs = 1:nruns;
end
SUBJECTS_DIR=getenv('SUBJECTS_DIR');
%% Set defaults
if ~exist('srcROI','var')
    srcROI = 'volume';
end
if ~exist('hemis','var')
    hemis = {'lh' 'rh'};
end
if ~exist('maps','var')
    %     maps = {'co' 'coecc' 'copol' 'copeakt' 'cosig1' 'cosig2' 'cosig3' 'cosig4' ...
    %         'var' 'varecc' 'varpol' 'varpeakt' 'var_as_co' 'varsig1' 'varsig2' 'varsig3' 'varsig4'};
    maps = {'co' 'coecc' 'copol' 'copeakt' 'cosig1' 'cosig2' 'cosig3' 'cosig4'};
end
disp(['Session_dir = ' session_dir]);
disp(['Found ' num2str(nruns) ' runs']);
disp(['Working on runs ' num2str(runs)]);
%% Average maps
progBar = ProgressBar(length(hemis),'Averaging maps...');
for hh = 1:length(hemis)
    clear sum_mat
    hemi = hemis{hh};
    for m = 1:length(maps)
        tmp = [];
        % Get the map values for all the runs
        for rr = runs
            tmpfile = fullfile(session_dir,'pRFs',d{rr},...
                [hemi '.' srcROI '.' maps{m} '.prfs']);
            if strcmp(srcROI,'cortex')
                lhtmpvol = load_nifti([tmpfile '.nii.gz']);
            else
                lhtmpvol = load_nifti([tmpfile '.orig.nii.gz']);
            end
            tmpvect = reshape(lhtmpvol.vol,...
                size(lhtmpvol.vol,1)*size(lhtmpvol.vol,2)*size(lhtmpvol.vol,3),1);
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
        save_nifti(nii,fullfile(session_dir,'pRFs',...
            [hemi '.' srcROI '.avg.' maps{m} '.prfs.nii.gz']));
        if strcmp(srcROI,'cortex') % Project to fsaverage_sym space
            sval = fullfile(session_dir,'pRFs',...
                [hemi '.' srcROI '.avg.' maps{m} '.prfs.nii.gz']);
            oval = fullfile(session_dir,'pRFs',...
                [hemi '.' srcROI '.avg.' maps{m} '.prfs.sym.nii.gz']);
            if strcmp(hemi,'lh')
                [~,~] = system(['mri_surf2surf --hemi lh --srcsubject ' subject_name ...
                    ' --sval ' sval ' --trgsubject fsaverage_sym --trgsurfval ' ...
                    oval]);
            elseif strcmp(hemi,'rh')
                [~,~] = system(['mri_surf2surf --hemi lh --srcsubject ' subject_name ...
                    '/xhemi --sval ' sval ' --trgsubject fsaverage_sym --trgsurfval ' ...
                    oval]);
                % average lh and rh
                lhvol = fullfile(session_dir,'pRFs',...
                    ['lh.' srcROI '.avg.' maps{m} '.prfs.sym.nii.gz']);
                rhvol = fullfile(session_dir,'pRFs',...
                    ['rh.' srcROI '.avg.' maps{m} '.prfs.sym.nii.gz']);
                lh = load_nifti(lhvol);
                rh = load_nifti(rhvol);
                tmp = [lh.vol,rh.vol];
                if strcmp('copol',maps{m}) || strcmp('varpol',maps{m})
                    tmp = [];
                    % Put right hemisphere prfs in left hemisphere coordinates
                    upper = rh.vol>0;
                    lower = rh.vol<0;
                    rh.vol(upper) = -rh.vol(upper) + pi;
                    rh.vol(lower) = -rh.vol(lower) - pi;
                    tmp = [lh.vol,rh.vol];
                    %tmpavg = nancirc_mean(tmp);
                    tmpavg = circ_mean(tmp');
                    tmpavg = tmpavg';
                    lh.vol = tmpavg;
                else
                    %tmpavg = nanmean(tmp,2);
                    tmpavg = mean(tmp,2);
                    lh.vol = tmpavg;
                end
                % Save nifti
                save_nifti(lh,fullfile(session_dir,'pRFs',...
                    ['mh.' srcROI '.avg.' maps{m} '.prfs.sym.nii.gz']));
                % Project back to subject space
                msval = fullfile(session_dir,'pRFs',...
                    ['mh.' srcROI '.avg.' maps{m} '.prfs.sym.nii.gz']);
                moval = fullfile(session_dir,'pRFs',...
                    ['mh.' srcROI '.avg.' maps{m} '.prfs.nii.gz']);
                [~,~] = system(['mri_surf2surf --hemi lh --srcsubject fsaverage_sym ' ...
                    ' --sval ' msval ' --trgsubject ' subject_name ' --trgsurfval ' ...
                    moval]);
            end
        else
            if strcmp(hemi,'rh') % need to complete lh first
                lh_vol_in = fullfile(session_dir,'pRFs',['lh.' srcROI '.avg.' maps{m} '.prfs.nii.gz']);
                rh_vol_in = fullfile(session_dir,'pRFs',['rh.' srcROI '.avg.' maps{m} '.prfs.nii.gz']);
                mh_vol_out = fullfile(session_dir,'pRFs',['mh.' srcROI '.avg.' maps{m} '.prfs.mni.nii.gz']);
                if strcmp('copol',maps{m}) || strcmp('varpol',maps{m})
                    rhtmpfile = fullfile(session_dir,'pRFs',['rh.' srcROI '.tmp.nii.gz']);
                    rh = load_nifti(rh_vol_in);
                    % Put right hemisphere prfs in left hemisphere coordinates
                    upper = rh.vol>0;
                    lower = rh.vol<0;
                    rh.vol(upper) = -rh.vol(upper) + pi;
                    rh.vol(lower) = -rh.vol(lower) - pi;
                    save_nifti(rh,rhtmpfile);
                    average_in_cvsMNI(session_dir,subject_name,maps{m},lh_vol_in,rhtmpfile,mh_vol_out)
                    delete(rhtmpfile);
                else
                    average_in_cvsMNI(session_dir,subject_name,maps{m},lh_vol_in,rh_vol_in,mh_vol_out)
                end
                % Project back to subject space
                out_vol = fullfile(session_dir,'pRFs',['mh.' srcROI '.avg.' maps{m} '.prfs.nii.gz']);
                ref_vol = fullfile(session_dir,'pRFs',['lh.' srcROI '.avg.' maps{m} '.prfs.nii.gz']);
                
                apply_cvs_inverse(mh_vol_out,out_vol,ref_vol,subject_name,SUBJECTS_DIR);
                
                %                 system(['mri_vol2vol --inv-morph --targ ' mh_vol_out ...
                %                     ' --m3z final_CVSmorph_tocvs_avg35_inMNI152.m3z' ...
                %                     ' --reg ' fullfile(SUBJECTS_DIR,subject_name,'cvs/subj2MNI.dat') ...
                %                     ' --mov ' ref_vol ' --o ' out_vol ' --inv-morph']);
            end
        end
    end
    progBar(hh);
end
disp('done.');
