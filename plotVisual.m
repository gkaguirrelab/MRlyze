function [outImage] = plotVisual(inVol,eccVol,polVol,roiInd,axLim,vArea)

% Plots values from an input volume in visual field coordinates
%
%   Usage:
%   [outImage] = plotVisual(inVol,eccVol,polVol,roiInd,axLim,vArea)
%
%   See also:
%   convert_image2surf
%
%   Written by Andrew S Bock Mar 2016

%% to do
% get the visual angle by eccentricity function for the pixels
%   combine ^ with the sigma by eccentricity function

%% set defaults
if ~exist('axLim','var')
    axLim = 90;
end
if ~exist('vArea','var')
    vArea = 'V1';
end
matSize = 200;
circPol = linspace(0,2*pi,matSize); % Polar angle sampling
cLines = linspace(0,log10(axLim),6); % Eccentricity lines
rLines = linspace(0,(2*pi) - (2*pi)/12,12); % Polar angle lines (spokes)
outImage = nan(matSize,matSize);
centerMat = (matSize+1)/2;
%% Load in volumes
% eccentricity data
if ischar(eccVol)
    if exist(eccVol,'file')
        [~,~,ext] = fileparts(eccVol);
        if strcmp(ext,'.gz')
            ecc = load_nifti(eccVol);
        elseif strcmp(ext,'.mgh') || strcmp(ext,'.mgz')
            ecc.vol = load_mgh(eccVol);
        end
    else
        error('data input does not correspond to an existing file');
    end
else
    ecc.vol = double(eccVol);
end
ecc = ecc.vol(roiInd);
% polar angle data
if ischar(polVol)
    if exist(polVol,'file')
        [~,~,ext] = fileparts(polVol);
        if strcmp(ext,'.gz')
            pol = load_nifti(polVol);
        elseif strcmp(ext,'.mgh') || strcmp(ext,'.mgz')
            pol.vol = load_mgh(polVol);
        end
    else
        error('data input does not correspond to an existing file');
    end
else
    pol.vol = double(polVol);
end
pol = pol.vol(roiInd);
% input data
if ischar(inVol)
    if exist(inVol,'file')
        [~,~,ext] = fileparts(inVol);
        if strcmp(ext,'.gz')
            in = load_nifti(inVol);
        elseif strcmp(ext,'.mgh') || strcmp(ext,'.mgz')
            in.vol = load_mgh(inVol);
        end
    else
        error('data input does not correspond to an existing file');
    end
else
    in.vol = double(inVol);
end
in = in.vol(roiInd);
%% Threshold based on axis limit
badind = ecc>axLim;
ecc(badind) = [];
pol(badind) = [];
in(badind) = [];
% in = icdf('normal',in); % convert p to z-score
% in(in==-inf) = -10;
% in(in==inf) = 10;
%% Convert polar to cartesian
[x,y] = pol2cart(pol,ecc);
% flip y
y = -y;
%% Make matrix
Sig = rf_ecc(ecc,vArea);
progBar = ProgressBar(matSize,'mixing paint...');
for i = 1:matSize % row (y)
    for j = 1:matSize % columns (x)
        tmpdist = sqrt( (i-centerMat)^2 + (j-centerMat)^2 );
        if tmpdist <= centerMat
            eccMat = 10^(tmpdist * log10(axLim)/(centerMat)); % log scale the value
            polMat = cart2pol(j-centerMat,i-centerMat);
            [matX,matY] = pol2cart(polMat,eccMat);
            tmpDist = sqrt( (matY - y).^2 + (matX - x).^2);
            tmpGauss = exp(-(tmpDist.^2)./(2*Sig.^2));
            tmpVals = tmpGauss.*in;
            outImage(i,j) = sum((tmpGauss .* tmpVals)) / sum(tmpGauss);
        end
    end
    if ~mod(i,10);progBar(i);end
end
%% Plot image
fullFigure;
subplot(1,1,1);%hold on;
axis off;
% bak = outImage;
%pImage = cdf('normal',outImage);
% threshImage = pImage;
% threshImage(threshImage>0.05) = nan;
% finalImage = log10(threshImage);
%finalImage = log10(pImage);
pcolor(outImage);
shading flat;
colormap(viridis);
axis off
axis square
hold on;
%% Create circles and spokes
clear cX cY
% Create circles
for i = 1:length(cLines)
    [cX(i,:),cY(i,:)] = pol2cart(circPol,(cLines(i)/max(cLines))*centerMat);
end
% Move to center, fix indexing
cX = cX + centerMat + 0.5;
cY = cY + centerMat + 0.5;
%% Plot circles and spokes
% Plot circles
for i = 1:length(cLines)
    plot(cX(i,:),cY(i,:),'k:','LineWidth',1);
    if cLines(i) ~= 0
        [tX,tY] = pol2cart(deg2rad(15),(cLines(i)/max(cLines))*(matSize-1)/2);
        text(-tX + max(cX(:))/2,tY + max(cY(:))/2,num2str(round(10^cLines(i))), 'FontSize', 20,...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end
% Plot spokes
for i = 1:length(rLines)
    [sX,sY] = pol2cart(rLines(i),(matSize-1)/2);
    plot([sX -sX] + max(cX(:))/2,[sY -sY] + max(cY(:))/2,'k:','LineWidth',1);
    [tX,tY] = pol2cart(rLines(i),(matSize-1)/2+5);
    text(tX + max(cX(:))/2,tY + max(cY(:))/2,num2str(rad2deg( rLines(i) )) , 'FontSize', 20,...
        'HorizontalAlignment','center','VerticalAlignment','middle');
end
axis square;
colorbar('NorthEastOutside');