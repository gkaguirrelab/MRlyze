function [varexp,B,BINT,R,RINT,STATS] = regress_fine_template(session_dir,runs,hemi,template,cluster)

% Run linear regression using a retinotopic template and movie correlation
% data
%
%   Usage:
%   [B,BINT,R,RINT,STATS] = regress_template(areas_template,movie_mat,corr_mat)
%
%   Written by Andrew S Bock Aug 2015

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
    corr_mat{i} = triu(corr_mat{i});
    movie_mat{i} = triu(movie_mat{i});
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
%% Fill in matrices
% Movie mat
movie_all = zero_mat;
movie_all(V1inds,V1inds) = movie_mat{1};
movie_all(V1inds,V2inds) = movie_mat{2};
movie_all(V1inds,V3inds) = movie_mat{3};
movie_all(V2inds,V2inds) = movie_mat{4};
movie_all(V2inds,V3inds) = movie_mat{5};
movie_all(V3inds,V3inds) = movie_mat{6};
% V1-V1
V{1} = zero_mat;
V{1}(V1inds,V1inds) = corr_mat{1} - mean(corr_mat{1}(:));
V{7} = zero_mat;
V{7}(V1inds,V1inds) = ones(size(corr_mat{1}));
V{7} = triu(V{7});
% V1-V2
V{2} = zero_mat;
V{2}(V1inds,V2inds) = corr_mat{2} - mean(corr_mat{2}(:));
V{8} = zero_mat;
V{8}(V1inds,V2inds) = ones(size(corr_mat{2}));
V{8} = triu(V{8});
% V1-V3
V{3} = zero_mat;
V{3}(V1inds,V3inds) = corr_mat{3} - mean(corr_mat{3}(:));
V{9} = zero_mat;
V{9}(V1inds,V3inds) = ones(size(corr_mat{3}));
V{9} = triu(V{9});
% V2-V2
V{4} = zero_mat;
V{4}(V2inds,V2inds) = corr_mat{4} - mean(corr_mat{4}(:));
V{10} = zero_mat;
V{10}(V2inds,V2inds) = ones(size(corr_mat{4}));
V{10} = triu(V{10});
% V2-V3
V{5} = zero_mat;
V{5}(V2inds,V3inds) =corr_mat{5} - mean(corr_mat{5}(:));
V{11} = zero_mat;
V{11}(V2inds,V3inds) = ones(size(corr_mat{5}));
V{11} = triu(V{11});
% V3-V3
V{6} = zero_mat;
V{6}(V3inds,V3inds) =corr_mat{6} - mean(corr_mat{6}(:));
V{12} = zero_mat;
V{12}(V3inds,V3inds) = ones(size(corr_mat{6}));
V{12} = triu(V{12});
%% Get upper triangular matrix and non-diagonals
diagind = zeros(size(movie_all));
diagind(1:length(movie_all)+1:end) = 1;
%goodind = diagind~=1 & movie_all~=0;
goodind = diagind~=1;
Y = movie_all(goodind);
tmpX = zeros(size(Y,1),12);
for i = 1:12
    tmp = V{i}(goodind);
    tmpX(:,i) = tmp - mean(tmp);
end
X = [ones(size(tmpX,1),1) tmpX];
%% Run regression
[B,BINT,R,RINT,STATS] = regress(Y,X);
%% Calc variance explained
for i = 1:12
    varexp(i) = var(B(i+1)*(X(:,i+1))) / var(Y); % since the first column is ones, use i+1 for X and B;
end
%% Save matrix
save(fullfile(pRF_dir,[template '.mat']),'Y','X','varexp','B','BINT','R',...
    'RINT','STATS','raw_corr_mat','raw_movie_mat','V1inds','V2inds','V3inds');
%% If running on cluster, save output value to text file
if cluster
    [~,~] = system(['echo ' num2str(varexp) ' > ' fullfile(pRF_dir,...
        [hemi '.' template '.varexp.txt'])]);
    [~,~] = system(['echo ' num2str(STATS(1)) ' > ' fullfile(pRF_dir,...
        [hemi '.' template '.rsquared.txt'])]);
    [~,~] = system(['echo ' num2str(B') ' > ' fullfile(pRF_dir,...
        [hemi '.' template '.betas.txt'])]);
end