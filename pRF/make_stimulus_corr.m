function [out_mat,srcind,trgind] = make_stimulus_corr(session_dir,runs,hemi,func,srcind,trgind,tcPart,leaveOut)

% Creates a cross-correlation matrix for stimulus driven data, averaging
%   across runs.
%
%   Usage:
%   [out_mat,srcind,trgind] = make_stimulus_corr(session_dir,runs,hemi,func,srcind,trgind,tcPart,leaveOut)
%
%   Written by Andrew S Bock Jul 2015

%% Find bold run directories
d = find_bold(session_dir);
if ~exist('tcPart','var')
    tcPart = 'full'; % 'H1' = 1st half; 'H2' = 2nd half
end
%% Read in occipital lable and functional volume
ct = 0;
for rr = 1:length(runs)
    run = runs(rr);
    func_vol = fullfile(session_dir,d{run},[func '.surf.0.1.inflated.' hemi '.nii.gz']);
    %func_vol = fullfile(session_dir,d{run},['s2.dbrf.tf.surf.0.1.inflated.' hemi '.nii.gz']);
    disp(['Loading ' func_vol '...']);
    tc = load_nifti(func_vol);
    disp('done.');
    tcs = squeeze(tc.vol)';
    switch tcPart
        case 'full'
            TRs = 1:size(tcs,1);
            srcTCs = tcs(TRs,srcind);
            trgTCs = tcs(TRs,trgind);
            tmp = fisher_z_corr(corr(srcTCs,trgTCs));
            tmp(tmp==inf) = 10;
            tmp(tmp==-inf) = -10;
            tmp_mat(rr,:,:) = tmp;
        case 'half'
            % 1st half
            ct = ct + 1;
            TRs1 = 1:(size(tcs,1)/2);
            if ~ismember(ct,leaveOut);
                srcTCs1 = tcs(TRs1,srcind);
                trgTCs1 = tcs(TRs1,trgind);
            else
                srcTCs1 = nan(length(TRs1),length(srcind));
                trgTCs1 = nan(length(TRs1),length(trgind));
            end
            tmp = fisher_z_corr(corr(srcTCs1,trgTCs1));
            tmp(tmp==inf) = 10;
            tmp(tmp==-inf) = -10;
            tmp_mat(ct,:,:) = tmp;
            % 2nd half
            ct = ct + 1;
            TRs2 = floor((1:(size(tcs,1)/2)) + size(tcs,1)/2);
            if ~ismember(ct,leaveOut);
                srcTCs2 = tcs(TRs2,srcind);
                trgTCs2 = tcs(TRs2,trgind);
            else
                srcTCs2 = nan(length(TRs2),length(srcind));
                trgTCs2 = nan(length(TRs2),length(trgind));
            end
            tmp = fisher_z_corr(corr(srcTCs2,trgTCs2));
            tmp(tmp==inf) = 10;
            tmp(tmp==-inf) = -10;
            tmp_mat(ct,:,:) = tmp;
    end
end
out_mat = squeeze(nanmean(tmp_mat,1));
disp('done.');