function submit_pRF(session_dir,rr,srcROI,hemi,srcind,datamat)

% Runs 'fit_pRF' to find the optimal pRF parameters, then save the output
% in a text file, in the following column format:
%
%   voxel index; pRF x-center; pRF y-center; pRF sigma; HRF time to peak
%
%   Usage:
%       submit_pRF(session_dir,rr,srcROI,hemi,srcind,datamat)
%
%   Written by Andrew S Bock May 2015

%% Find bold run directories
d = listdir(fullfile(session_dir,'*BOLD_*'),'dirs');
if isempty(d)
    d = listdir(fullfile(session_dir,'*bold_*'),'dirs');
end
if isempty(d)
    d = listdir(fullfile(session_dir,'*EPI_*'),'dirs');
end
if isempty(d)
    d = listdir(fullfile(session_dir,'RUN*'),'dirs');
end
%% get the source timecourse
if strcmp(srcROI,'occipital') || strcmp(srcROI,'cortex')
    srcfile = fullfile(session_dir,d{rr},['sdbrf.tf_surf.' hemi '.nii.gz']);
elseif strcmp(srcROI,'volume')
    srcfile = fullfile(session_dir,d{rr},'sdbrf.tf.nii.gz');
end
src = load_nifti(srcfile);
srcdims = size(src.vol);
srctc = reshape(src.vol,srcdims(1)*srcdims(2)*srcdims(3),srcdims(4))';
srctc = srctc(:,srcind);
%% Get TR, X, Y, images, init_params, lower_bound, upper_bound
load(datamat);
tic
[opt_params,fval] = fit_pRF(srctc,TR,X,Y,images,init_params(srcind,:));
toc
%% Save opt_params in a output text file
if strcmp(srcROI,'occipital') || strcmp(srcROI,'cortex')
    out_file = fullfile(session_dir,['Run_' num2str(rr) '_' srcROI '_' hemi '_opt_params.txt']);
elseif strcmp(srcROI,'volume')
    out_file = fullfile(session_dir,['Run_' num2str(rr) '_' srcROI '_opt_params.txt']);
end
f = fopen(out_file,'a');
fprintf(f, [num2str([srcind opt_params fval]) '\n']);