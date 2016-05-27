function plotVisual(inVol,eccVol,polVol,roiInd,axLim,vArea)

% Plots values from an input volume in visual field coordinates
%
%   Usage:
%   plotVisual(inVol,eccVol,polVol,sigVol,roiInd,roiName,axLim)
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
pscale = 100; % pixel scaling factor
circPol = 0:pi/50:2*pi; % Polar angle sampling
cLines = linspace(0,log10(axLim),6); % Eccentricity lines
rLines = linspace(0,(2*pi) - (2*pi)/12,12); % Polar angle lines (spokes)
pLabels = pscale*(cLines(end) + (cLines(end) - cLines(end-1))/4); % Polar angle label eccentricity
%% Set display
fullFigure;
subplot(1,1,1);%hold on;
axis off;
%% Load in volumes
ecc = load_nifti(eccVol);
ecc = ecc.vol(roiInd);
pol = load_nifti(polVol);
pol = pol.vol(roiInd);
in = load_nifti(inVol);
in = in.vol(roiInd);
%% Remove values outside axLim
badind = ecc>axLim;
ecc(badind) = [];
pol(badind) = [];
in(badind) = [];
%% Log scale ecc
logecc = log10(ecc)*pscale;
logecc(logecc<0) = 0;
%% Convert polar to cartesian
[x,y] = pol2cart(pol,logecc);
% flip y
y = -y;
%% Create circles and spokes
% Create circles
for i = 1:length(cLines)
    [cX(i,:),cY(i,:)] = pol2cart(circPol,cLines(i)*pscale);
end
% Make all postivie
cX = cX - min(cX(:));
cY = cY - min(cY(:));
% Get circle x, y
xCenter = max(cX(:))/2;
yCenter = max(cY(:))/2;
% Get pixel size by degrees eccentricity
maxPix = max(cX(:))/2;
%% Plot Background
patch('XData',cX(end,:),'YData',cY(end,:),...
    'FaceColor',[1 1 1],'HandleVisibility', 'off');
hold on;
%% Create image
tmpImage = zeros(ceil(max(cX(:))),ceil(max(cY(:))));
tmpX = round(x + max(cX(:))/2);
tmpY = round(y + max(cY(:))/2);
xMat = repmat(1:size(tmpImage,1),size(tmpImage,2),1);
yMat = repmat(1:size(tmpImage,2),size(tmpImage,1),1)';
valMat = nan([size(tmpImage,1)*size(tmpImage,2),length(x)]);
weightMat = nan([size(tmpImage,1)*size(tmpImage,2),length(x)]);
% Make the matrix of values and weights
progBar = ProgressBar(length(tmpX)','Making image matrices...');
for i = 1:length(tmpX)
    eccSig = rf_ecc(ecc(i),vArea);
    pixSig = (log10(eccSig))*pscale;
    if pixSig < .1
        pixSig = .1;
    end   
    tmpDist = sqrt( (yMat - tmpY(i)).^2 + (xMat - tmpX(i)).^2);
    tmpGauss = exp(-(tmpDist(:).^2)/(2*pixSig.^2));
    valMat(:,i) = tmpGauss*in(i);
    weightMat(:,i) = tmpGauss;
    if ~mod(i,100);progBar(i);end
end
%% Get weighted average
outImage = nan(size(weightMat,1),1);
progBar = ProgressBar(length(weightMat),'Calculating weighted average...');
for i = 1:length(weightMat)
    tmpVals = valMat(i,:);
    tmpWeights = weightMat(i,:);
    outImage(i) = sum((tmpWeights .* tmpVals)) / sum(tmpWeights);
    if ~mod(i,100);progBar(i);end
end
%% Make Final Image
finalImage = reshape(outImage,size(tmpImage));
for i = 1:size(finalImage,1)
    for j = 1:size(finalImage,2)
        matDist = sqrt((i-xCenter)^2 + (j-yCenter)^2);
        if matDist >= xCenter-1
            finalImage(i,j) = nan;
        end
    end
end
finalImage = log10(finalImage);
%% Plot image
pcolor(finalImage);
shading flat;
colormap(flipud(hot(2000)));
hold on;
%% Plot polar grid
% Plot circles
for i = 1:length(cLines)
    plot(cX(i,:),cY(i,:),'k:','LineWidth',1);
    if cLines(i) ~= 0
        [tX,tY] = pol2cart(deg2rad(15),cLines(i)*pscale);
        text(-tX + max(cX(:))/2,tY + max(cY(:))/2,num2str(round(10^cLines(i))), 'FontSize', 20,...
            'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor','white');
    end
end
% Plot spokes
for i = 1:length(rLines)
    [sX,sY] = pol2cart(rLines(i),cLines(end)*pscale);
    plot([sX -sX] + max(cX(:))/2,[sY -sY] + max(cY(:))/2,'k:','LineWidth',1);
    [tX,tY] = pol2cart(rLines(i),pLabels);
    text(tX + max(cX(:))/2,tY + max(cY(:))/2,num2str(rad2deg( rLines(i) )) , 'FontSize', 20,...
        'HorizontalAlignment','center','VerticalAlignment','middle');
    %     text(-tX,-tY,num2str(rad2deg( mod(pi+rLines(i),2*pi) )) ,'FontSize', 20,...
    %         'HorizontalAlignment','center','VerticalAlignment','middle');
end
axis square;
colorbar('EastOutside');