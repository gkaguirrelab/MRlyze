function [R2,B] = regress_template_V2V3(session_dir,saveDir,templateType,runs,hemi,func,template,cluster,tcPart,leaveOut)

% Run linear regression using a retinotopic template and movie correlation
% data
%
%   Usage:
%   [R2,B] = regress_template_V2V3(session_dir,saveDir,templateType,runs,hemi,func,template,cluster,tcPart,leaveOut)
%
%   Written by Andrew S Bock Aug 2015

%% Set defaults
areas = {'V1' 'V2' 'V3'};
switch templateType
    case 'coarse'
        pRF_dir = fullfile(session_dir,'pRFs','coarse_model_templates','decimated_templates');
    case 'fine'
        pRF_dir = fullfile(session_dir,'pRFs','fine_model_templates','decimated_templates');
    case 'pRF'
        pRF_dir = fullfile(session_dir,'pRFs','pRF_templates','decimated_templates');
    case 'anat'
        pRF_dir = fullfile(session_dir,'anat_templates','decimated_templates');
end
if ~exist('tcPart','var')
    tcPart = 'full'; % 'H1' = 1st half; 'H2' = 2nd half
end
%% Get template and movie correlation matrices
ct = 0;
for ss = 1:length(areas)
    for tt = ss:length(areas)
        ct = ct + 1;
        areas_srctemplate = fullfile(pRF_dir,[hemi '.areas.' template '.nii.gz']);
        ecc_srctemplate = fullfile(pRF_dir,[hemi '.ecc.' template '.nii.gz']);
        pol_srctemplate = fullfile(pRF_dir,[hemi '.pol.' template '.nii.gz']);
        [corr_mat{ct},srcind,trgind] = ...
            make_model_corr(areas{ss},areas{tt},areas_srctemplate,ecc_srctemplate,pol_srctemplate);
        [movie_mat{ct}] = make_stimulus_corr(session_dir,runs,hemi,func,srcind,trgind,tcPart,leaveOut);
    end
end
%% Create zero full matrix
V1length = size(corr_mat{1},2);
V2length = size(corr_mat{2},2);
V3length = size(corr_mat{3},2);
rowV1inds = 1:V1length;
rowV2inds = 1+V1length:V1length+V2length;
colV2inds = 1:V2length;
colV3inds = 1+V2length:V2length+V3length;
zero_mat = zeros(V1length+V2length,V2length+V3length);
zero_mat(rowV2inds,colV2inds) = nan; % 'nan' use to remove V2-V2 later
%% Regress
% Zero mats
corrV1V2 = zero_mat;
mcorrV1V2 = zero_mat;
corrV1V3 = zero_mat;
mcorrV1V3 = zero_mat;
corrV2V3 = zero_mat;
mcorrV2V3 = zero_mat;
movieMat = zero_mat;
% V1-V2
corrV1V2(rowV1inds,colV2inds) = corr_mat{2};
corrV1V2(isnan(corrV1V2)) = []; % remove V2-V2
corrV1V2 = corrV1V2(:) - mean(corrV1V2(:));
mcorrV1V2(rowV1inds,colV2inds) = 1;
mcorrV1V2(isnan(mcorrV1V2)) = []; % remove V2-V2
mcorrV1V2 = mcorrV1V2(:);
movieMat(rowV1inds,colV2inds) = movie_mat{2};
% V1-V3
corrV1V3(rowV1inds,colV3inds) = corr_mat{3};
corrV1V3(isnan(corrV1V3)) = [];
corrV1V3 = corrV1V3(:) - mean(corrV1V3(:));
mcorrV1V3(rowV1inds,colV3inds) = 1;
mcorrV1V3(isnan(mcorrV1V3)) = [];
mcorrV1V3 = mcorrV1V3(:);
movieMat(rowV1inds,colV3inds) = movie_mat{3};
% V2-V3
corrV2V3(rowV2inds,colV3inds) = corr_mat{5};
corrV2V3(isnan(corrV2V3)) = [];
corrV2V3 = corrV2V3(:) - mean(corrV2V3(:));
mcorrV2V3(rowV2inds,colV3inds) = 1;
mcorrV2V3(isnan(mcorrV2V3)) = [];
mcorrV2V3 = mcorrV2V3(:);
movieMat(rowV2inds,colV3inds)= movie_mat{5};
% Remove V2-V2
movieMat(isnan(movieMat)) = [];
% Regression
X = [corrV1V2,corrV1V3,corrV2V3,mcorrV1V2,mcorrV1V3,mcorrV2V3];
Y = movieMat(:);
[B,~,~,~,STATS] = regress(Y,X);
R2 = STATS(1);
for i = 1:length(B)
    varexp(i) = var(B(i)*X(:,i)) / var(Y);
end
%% Save text files
if cluster
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
    [~,~] = system(['echo ' num2str(varexp) ' > ' fullfile(saveDir,...
        [hemi '.' template '.varexp.txt'])]);
    [~,~] = system(['echo ' num2str(R2) ' > ' fullfile(saveDir,...
        [hemi '.' template '.rsquared.txt'])]);
    [~,~] = system(['echo ' num2str(B') ' > ' fullfile(saveDir,...
        [hemi '.' template '.betas.txt'])]);
end