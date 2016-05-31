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
    axLim = 37;
end
if ~exist('vArea','var')
    vArea = 'V1';
end
weightThresh = 1;
matSize = 501;
circPol = linspace(0,2*pi,matSize); % Polar angle sampling
cLines = linspace(0,log10(axLim),6); % Eccentricity lines
rLines = linspace(0,(2*pi) - (2*pi)/12,12); % Polar angle lines (spokes)
pLabels = matSize*(cLines(end) + (cLines(end) - cLines(end-1))/4); % Polar angle label eccentricity
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
%% Threshold based on axis limit
badind = ecc>axLim;
ecc(badind) = [];
pol(badind) = [];
in(badind) = [];
in = icdf('normal',in); % convert p to z-score
in(in==-inf) = -10;
in(in==inf) = 10;
%% Convert polar to cartesian
[x,y] = pol2cart(pol,ecc);
% flip y
y = -y;
%% Make matrix
outImage = nan(matSize,matSize);
centerMat = [round(matSize/2) round(matSize/2)];
Sig = rf_ecc(ecc,vArea);
progBar = ProgressBar(matSize,'mixing paint...');
for i = 1:matSize % row (y)
    for j = 1:matSize % columns (x)
        tmpdist = sqrt( (i-centerMat(1))^2 + (j-centerMat(2))^2 );
        if tmpdist <= centerMat(1)-1
            eccMat = 10^(tmpdist * log10(axLim)/(centerMat(1)-1)); % log scale the value
            polMat = cart2pol(j-centerMat(2),i-centerMat(1));
            [matX,matY] = pol2cart(polMat,eccMat);
            tmpDist = sqrt( (matY - y).^2 + (matX - x).^2);
            tmpGauss = exp(-(tmpDist.^2)./(2*Sig.^2));
            tmpGauss(tmpGauss<0.5) = 0;
            tmpVals = tmpGauss.*in;
            if sum(tmpGauss) >= weightThresh
                outImage(i,j) = sum((tmpGauss .* tmpVals)) / sum(tmpGauss);
            end
        end
    end
    progBar(i);
end
%% Plot image
% bak = outImage;
pImage = cdf('normal',outImage);
% threshImage = pImage;
% threshImage(threshImage>0.05) = nan;
% finalImage = log10(threshImage);
finalImage = log10(pImage);
pcolor(finalImage);
shading flat;
colormap(flipud(hot(2000)));
axis off
axis square
caxis([-5 -0]);
hold on;
%% Create circles and spokes
% Create circles
for i = 1:length(cLines)
    [cX(i,:),cY(i,:)] = pol2cart(circPol,cLines(i)*matSize/2);
end
% Make all postivie
cX = cX - min(cX(:));
cY = cY - min(cY(:));
% Get circle x, y
xCenter = max(cX(:))/2;
yCenter = max(cY(:))/2;
%% Plot Background
patch('XData',cX(end,:),'YData',cY(end,:),...
    'FaceColor',[1 1 1],'HandleVisibility', 'off');
hold on;








%%




%%


%%



%%



%%
% scale and move x and y
x = round(x) + round(size(imageMat,2)/2);
y = round(y) + round(size(imageMat,1)/2);
ct = 0;
for i = 1:length(in)
    if isnan(imageMat(y(i),x(i)))
    imageMat(y(i),x(i)) = in(i);
    else
        ct = ct + 1;
    end
end




%%
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
    if eccSig < 1
        eccSig = 1;
    end
    pixSig = (log10(eccSig))*pscale;
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