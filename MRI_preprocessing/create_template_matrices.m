function [obs_all,pred_all,srcind,trgind,rowV1inds,rowV2inds,rowV3inds,colV1inds,colV2inds,colV3inds] = create_template_matrices(session_dir,pRF_dir,template,hemi,func,runs,templateSize)

% Creates the observed and predicted the cross-correlation matrices which
%   result from the retinotopic template fitting pipeline.
%
%   Written by Andrew S Bock Jan 2016

%% Set defaults
if ~exist('hemi','var')
    hemi = 'lh';
end
if ~exist('func','var')
    func = 's5.dbrf.tf';
end
if ~exist('runs','var')
    runs = [2 4 6];
end
if ~exist('templateSize','var')
    templateSize = 'full';
end
areas = {'V1' 'V2' 'V3'};
%% Define templates
areas_srctemplate = fullfile(pRF_dir,[hemi '.areas.' template '.nii.gz']);
ecc_srctemplate = fullfile(pRF_dir,[hemi '.ecc.' template '.nii.gz']);
pol_srctemplate = fullfile(pRF_dir,[hemi '.pol.' template '.nii.gz']);
% Get template and movie correlation matrices
ct = 0;
for ii = 1:length(areas)
    for jj = ii:length(areas)
        ct = ct + 1;
        [pred_mat{ct},srcind{ct},trgind{ct}] = ...
            make_model_corr(areas{ii},areas{jj},areas_srctemplate,ecc_srctemplate,pol_srctemplate);
        [obs_mat{ct}] = make_stimulus_corr(session_dir,runs,hemi,func,srcind{ct},trgind{ct});
    end
end
% Create Movie full matrix
V1length = size(obs_mat{1},2);
V2length = size(obs_mat{2},2);
V3length = size(obs_mat{3},2);
% Fill in matrices
switch templateSize
    case 'full'
        rowV1inds = 1:V1length;
        rowV2inds = 1+V1length:V1length+V2length;
        rowV3inds = 1+V1length+V2length:V1length+V2length+V3length;
        colV1inds = 1:V1length;
        colV2inds = 1+V1length:V1length+V2length;
        colV3inds = 1+V1length+V2length:V1length+V2length+V3length;
        full_length = V1length + V2length + V3length;
        zero_mat = zeros(full_length);
        % Movie mat
        obs_all = zero_mat;
        obs_all(rowV1inds,colV1inds) = obs_mat{1};
        obs_all(rowV1inds,colV2inds) = obs_mat{2};
        obs_all(rowV2inds,colV1inds) = obs_mat{2}';
        obs_all(rowV1inds,colV3inds) = obs_mat{3};
        obs_all(rowV3inds,colV1inds) = obs_mat{3}';
        obs_all(rowV2inds,colV2inds) = obs_mat{4};
        obs_all(rowV2inds,colV3inds) = obs_mat{5};
        obs_all(rowV3inds,colV2inds) = obs_mat{5}';
        obs_all(rowV3inds,colV3inds) = obs_mat{6};
        % nan out diagonal
        for i = 1:length(obs_all);obs_all(i,i) = nan;end
    case 'V2_V3'
        rowV1inds = 1:V1length;
        rowV2inds = 1+V1length:V1length+V2length;
        rowV3inds = [];
        colV1inds = [];
        colV2inds = 1:V2length;
        colV3inds = 1+V2length:V2length+V3length;
        %zero_mat = zeros(V1length,V2length+V3length);
        zero_mat = zeros(V1length+V2length,V2length+V3length);
        % Movie mat
        obs_all = zero_mat;
        obs_all(rowV1inds,colV2inds) = obs_mat{2};
        obs_all(rowV1inds,colV3inds) = obs_mat{3};
        obs_all(rowV2inds,colV3inds) = obs_mat{5};
        
end
% Create template full matrix
V1length = size(pred_mat{1},2);
V2length = size(pred_mat{2},2);
V3length = size(pred_mat{3},2);
% Fill in matrices
switch templateSize
    case 'full'
        rowV1inds = 1:V1length;
        rowV2inds = 1+V1length:V1length+V2length;
        rowV3inds = 1+V1length+V2length:V1length+V2length+V3length;
        colV1inds = 1:V1length;
        colV2inds = 1+V1length:V1length+V2length;
        colV3inds = 1+V1length+V2length:V1length+V2length+V3length;
        full_length = V1length + V2length + V3length;
        zero_mat = zeros(full_length);
        % Movie mat
        pred_all = zero_mat;
        pred_all(rowV1inds,colV1inds) = pred_mat{1};
        pred_all(rowV1inds,colV2inds) = pred_mat{2};
        pred_all(rowV2inds,colV1inds) = pred_mat{2}';
        pred_all(rowV1inds,colV3inds) = pred_mat{3};
        pred_all(rowV3inds,colV1inds) = pred_mat{3}';
        pred_all(rowV2inds,colV2inds) = pred_mat{4};
        pred_all(rowV2inds,colV3inds) = pred_mat{5};
        pred_all(rowV3inds,colV2inds) = pred_mat{5}';
        pred_all(rowV3inds,colV3inds) = pred_mat{6};
        % nan out diagonal
        for i = 1:length(pred_all);pred_all(i,i) = nan;end
    case 'V2_V3'
        rowV1inds = 1:V1length;
        rowV2inds = 1+V1length:V1length+V2length;
        rowV3inds = [];
        colV1inds = [];
        colV2inds = 1:V2length;
        colV3inds = 1+V2length:V2length+V3length;
        %zero_mat = zeros(V1length,V2length+V3length);
        zero_mat = zeros(V1length+V2length,V2length+V3length);
        % Movie mat
        pred_all = zero_mat;
        pred_all(rowV1inds,colV2inds) = pred_mat{2};
        pred_all(rowV1inds,colV3inds) = pred_mat{3};
        pred_all(rowV2inds,colV3inds) = pred_mat{5};
end