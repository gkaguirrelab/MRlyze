function [corr_mat,srcind,trgind,srcecc,trgecc] = make_model_corr(src_area,trg_area,areas_srctemplate,ecc_srctemplate,pol_srctemplate)

% Creates a cross-correlation matrix for retinotopic models
%
%   Usage:
%   [corr_mat,srcind,trgind,srcecc,trgecc] = make_model_corr(src_area,trg_area,areas_srctemplate,ecc_srctemplate,pol_srctemplate)
%
%   Written by Andrew S Bock Jul 2015

%% set defaults
if ~exist('src_area','var')
    src_area = 'V1';
end
if ~exist('trg_area','var')
    trg_area = 'V1';
end
%% Read in pRF area file, get ROI indices
areas = load_nifti(areas_srctemplate);
ECC = load_nifti(ecc_srctemplate);
POL = load_nifti(pol_srctemplate);
% source indicies
if strcmp(src_area,'V1')
    tmpi = find(areas.vol == -1 | areas.vol == 1);
elseif strcmp(src_area,'V2')
    tmpi = find(areas.vol == -2 | areas.vol == 2);
elseif strcmp(src_area,'V3')
    tmpi = find(areas.vol == -3 | areas.vol == 3);
elseif strcmp(src_area,'V1V2')
    tmpi = find(areas.vol >= -2 & areas.vol <= 2);
elseif strcmp(src_area,'Vall')
    tmpi = find(areas.vol >= -3 & areas.vol <= 3);
end
[~,tmpsi] = sort(ECC.vol(tmpi));
srcind = tmpi(tmpsi);
% target indicies
if strcmp(trg_area,'V1')
    tmpi = find(areas.vol == -1 | areas.vol == 1);
elseif strcmp(trg_area,'V2')
    tmpi = find(areas.vol == -2 | areas.vol == 2);
elseif strcmp(trg_area,'V3')
    tmpi = find(areas.vol == -3 | areas.vol == 3);
elseif strcmp(trg_area,'V1V2')
    tmpi = find(areas.vol >= -2 & areas.vol <= 2);
elseif strcmp(trg_area,'Vall')
    tmpi = find(areas.vol >= -3 & areas.vol <= 3);
end
[~,tmpsi] = sort(ECC.vol(tmpi));
trgind = tmpi(tmpsi);
%% Get Ecc and Pol
srcecc = ECC.vol(srcind);
srcpol = POL.vol(srcind);
trgecc = ECC.vol(trgind);
trgpol = POL.vol(trgind);
[srcx,srcy] = pol2cart(srcpol,srcecc);
[trgx,trgy] = pol2cart(trgpol,trgecc);
allind = [srcind' trgind'];
allecc = [srcecc' trgecc'];
allx = [srcx' trgx']';
ally = [srcy' trgy']';
%% Make stimulus movie
% x is going to be the distance (x,y) to the other points
numTRs = 1000;
randMov = randn(numTRs,size(allecc,2));
stimMov = randMov;
%stimMov = nan(numTRs,size(allecc,2));
vertMov = nan(numTRs,size(allecc,2));
disp(['making fake data ' src_area ' to ' trg_area '...']);
% for yy = 1:size(allind,2)
%     yydists = sqrt( (allx(yy) - allx).^2 + (ally(yy) - ally).^2 ); % in visual angle
%     % Create kernel based on stimulus correlation
%     P = [-0.000254159842224    0.008377945550527   -0.099740350218546    0.843457690087762];% based on actual Movie data, in pixels
%     stimKernel = polyval(P,yydists); % ~10 pixels/ visual angle for the movie data kernel code
%     stimKernel(stimKernel<0) = 0;
%     %stimKernel = exp(1).^(-yydists/2);
%     %     stimKernel(stimKernel == Inf) = nan; % Set Inf values to non-Inf max values
%     %     stimKernel(isnan(stimKernel)) = max(stimKernel);
%     %     stimKernel = ((stimKernel - min(stimKernel)) / max(stimKernel - min(stimKernel)));
%     % Multiply by receptive field location kernel
%     stimMov(:,yy) = randMov*stimKernel;
% end
progBar = ProgressBar(size(allind,2),'looping through vertices...');
for xx = 1:size(allind,2)
    xxdists = sqrt((allx(xx) - allx).^2 + (ally(xx) - ally).^2);
    % Create kernel based on receptive field locations
    areaind = areas.vol(allind(xx));
    if areaind == -1 || areaind == 1
        vtx_area = 'V1';
    elseif areaind == -2 || areaind == 2
        vtx_area = 'V2';
    elseif areaind == -3 || areaind == 3
        vtx_area = 'V2';
    end
    center_sig = rf_ecc(allecc(xx),vtx_area);
    center_rfKernel = exp(-(xxdists.^2)/(2*center_sig.^2));
    rfKernel = center_rfKernel;
    % add a Gaussian surround (DoG model)
    %     surround_sig = rf_ecc_surround(allecc(xx),vtx_area);
    %     surround_index = rf_ecc_surround_index(allecc(xx),vtx_area);
    %     surround_amp = surround_index ./ (surround_sig ./ center_sig);
    %     surround_amp(surround_amp<0) = 0;
    %     surround_rfKernel = exp(-(xxdists.^2)/(2*surround_sig.^2));
    %     rfKernel = center_rfKernel - surround_rfKernel*surround_amp;
    vertMov(:,xx) = stimMov*rfKernel;
    if ~mod(xx,1000);progBar(xx);end
end
disp('done.');
%% Make the correlation matrix (V1)
tmp = fisher_z_corr(corr(vertMov(:,1:length(srcind)),vertMov(:,length(srcind)+1:end)));
tmp(tmp==inf) = 10;
tmp(tmp==-inf) = -10;
%corr_mat = tmp;
%corr_mat = tmp-.2; % negative 0.2 offset
%corr_mat = corr(vertMov(:,1:length(srcind)),vertMov(:,length(srcind)+1:end));