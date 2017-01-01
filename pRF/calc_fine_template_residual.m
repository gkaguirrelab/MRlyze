function [residualSS,V1inds,V2inds,V3inds] = calc_fine_template_residual(session_dir,runs,hemi,template)

% Calculates the residual sum of squares between a template correlation 
%   matrix and a movie correlation matrix on the decimated surface.
%
%   Usage:
%   [residualSS,V1inds,V2inds,V3inds] = calc_fine_template_residual(session_dir,runs,hemi,template)
%
%   Written by Andrew S Bock Oct 2015

%% Set defaults
areas = {'V1' 'V2' 'V3'};
pRF_dir = fullfile(session_dir,'pRFs','fine_model_templates','decimated_templates');

%% Get template and movie correlation matrices
ct = 0;
for ss = 1:length(areas)
    for tt = ss:length(areas)
        ct = ct + 1;
        areas_srctemplate = fullfile(pRF_dir,[hemi '.areas.' template '.nii.gz']);
        ecc_srctemplate = fullfile(pRF_dir,[hemi '.ecc.' template '.nii.gz']);
        pol_srctemplate = fullfile(pRF_dir,[hemi '.pol.' template '.nii.gz']);
        [raw_corr_mat{ct},srcind,trgind] = ...
            make_model_corr(areas{ss},areas{tt},areas_srctemplate,ecc_srctemplate,pol_srctemplate);
        [raw_movie_mat{ct}] = make_stimulus_corr(session_dir,runs,hemi,srcind,trgind);
    end
end
%% Extract upper triangular part of matrices
corr_mat = raw_corr_mat;
movie_mat = raw_movie_mat;
for i = [1,4,6]
    corr_mat{i} = triu(raw_corr_mat{i});
    movie_mat{i} = triu(raw_movie_mat{i});
end
%% Create zero full matrix
V1length = size(corr_mat{1},2);
V2length = size(corr_mat{2},2);
V3length = size(corr_mat{3},2);
V1inds = 1:V1length;
V2inds = 1+V1length:V1length+V2length;
V3inds = 1+V1length+V2length:V1length+V2length+V3length;
full_length = V1length + V2length + V3length;
zero_mat = zeros(full_length);
%% Create full matrices
% Movie matrices
raw_movie_all = zero_mat;
raw_movie_all(V1inds,V1inds) = raw_movie_mat{1};
raw_movie_all(V1inds,V2inds) = raw_movie_mat{2};
raw_movie_all(V2inds,V1inds) = raw_movie_mat{2}';
raw_movie_all(V1inds,V3inds) = raw_movie_mat{3};
raw_movie_all(V3inds,V1inds) = raw_movie_mat{3}';
raw_movie_all(V2inds,V2inds) = raw_movie_mat{4};
raw_movie_all(V2inds,V3inds) = raw_movie_mat{5};
raw_movie_all(V3inds,V2inds) = raw_movie_mat{5}';
raw_movie_all(V3inds,V3inds) = raw_movie_mat{6};
% Correlation matrices
raw_corr_all = zero_mat;
raw_corr_all(V1inds,V1inds) = raw_corr_mat{1};
raw_corr_all(V1inds,V2inds) = raw_corr_mat{2};
raw_corr_all(V2inds,V1inds) = raw_corr_mat{2}';
raw_corr_all(V1inds,V3inds) = raw_corr_mat{3};
raw_corr_all(V3inds,V1inds) = raw_corr_mat{3}';
raw_corr_all(V2inds,V2inds) = raw_corr_mat{4};
raw_corr_all(V2inds,V3inds) = raw_corr_mat{5};
raw_corr_all(V3inds,V2inds) = raw_corr_mat{5}';
raw_corr_all(V3inds,V3inds) = raw_corr_mat{6};
%% Loop through all vertices
residualSS = nan(4,full_length);
for i = 1:full_length
    tmpMov = raw_movie_all(:,i);
    tmpCorr = raw_corr_all(:,i);
    tmpMov(i) = [];
    tmpCorr(i) = [];
    tmpV1mov = raw_movie_all(V1inds,i);
    tmpV1corr = raw_corr_all(V1inds,i);
    if ismember(i,V1inds)
        j = find(ismember(V1inds,i));
        tmpV1mov(j) = [];
        tmpV1corr(j) = [];
    end
    tmpV2mov = raw_movie_all(V2inds,i);
    tmpV2corr = raw_corr_all(V2inds,i);
    if ismember(i,V2inds)
        j = find(ismember(V2inds,i));
        tmpV2mov(j) = [];
        tmpV2corr(j) = [];
    end
    tmpV3mov = raw_movie_all(V3inds,i);
    tmpV3corr = raw_corr_all(V3inds,i);
    if ismember(i,V3inds)
        j = find(ismember(V3inds,i));
        tmpV3mov(j) = [];
        tmpV3corr(j) = [];
    end
    residualSS(1,i) = sum( (tmpV1mov - tmpV1corr) .^2);
    residualSS(2,i) = sum( (tmpV2mov -tmpV2corr) .^2);
    residualSS(3,i) = sum( (tmpV3mov -tmpV3corr) .^2);
    residualSS(4,i) = sum( (tmpMov - tmpCorr) .^2);
end
%% Save as .mat file
save(fullfile(pRF_dir,[hemi '.' template '.residualSS.mat']),'residualSS','V1inds','V2inds','V3inds');