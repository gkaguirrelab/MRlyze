function plotVisual(inVol,eccVol,polVol,sigVol,roiInd,roiName,axLim,caxLims)

% Plots values from an input volume in visual field coordinates
%
%   Usage:
%   plotVisual(inVol,eccVol,polVol,sigVol,roiInd,roiName,axLim)
%
%   See also:
%   convert_image2surf
%
%   Written by Andrew S Bock Mar 2016
%% set defaults
if ~exist('roiName','var') || isempty(roiName)
    roiName = 'V1';
end
if ~exist('axLim','var')
    axLim = 90;
end
if ~exist('caxLims','var')
    caxLims = [0 0.05];
end
pscale = 100;
circPol = 0:pi/50:2*pi; % Polar angle sampling
circEcc = ones(size(circPol)); % Eccentricity sampling
cLines = linspace(0,log10(axLim),6); % Eccentricity lines
rLines = linspace(0,(2*pi) - (2*pi)/12,12); % Polar angle lines (spokes)
pLabels = pscale*(cLines(end) + (cLines(end) - cLines(end-1))/2); % Polar angle label eccentricity
%% Set display
[~,screensize] = fullFigure;
subplot(1,1,1);%hold on;
axis off;
display.distance = 63; % assume 63 cm (average arm length)
display.width = 32.5; % assume 32.5 cm (laptop width)
display.resolution = screensize([3 4]);
%% Load in volumes
ecc = load_nifti(eccVol);
ecc = ecc.vol(roiInd);
logecc = log10(ecc);
logecc(logecc<0) = 0;
pol = load_nifti(polVol);
pol = pol.vol(roiInd);
if isempty(sigVol)
    sig = rf_ecc(ecc,roiName);
else
    sig = load_nifti(sigVol);
    sig = sig.vol(roiInd);
end
% load input volume
in = load_nifti(inVol);
in = in.vol(roiInd);
%% Convert polar to cartesian
[x,y] = pol2cart(pol,logecc*pscale);
% flip y
y = -y;
%% Plot circles and spokes
% Create circles
for i = 1:length(cLines)
    [cX(i,:),cY(i,:)] = pol2cart(circPol,cLines(i)*pscale);
end
% Make all postivie
cX = cX - min(cX(:));
cY = cY - min(cY(:));
% White background
patch('XData',cX(end,:),'YData',cY(end,:),...
    'FaceColor',[1 1 1],'HandleVisibility', 'off');
hold on;
tmpImage = zeros(ceil(max(cX(:)) - min(cX(:))),ceil(max(cY(:)) - min(cY(:))));
tmpX = round(x + max(cX(:))/2);
tmpY = round(y + max(cY(:))/2);
xMat = repmat(1:size(tmpImage,1),size(tmpImage,2),1);
yMat = repmat(1:size(tmpImage,2),size(tmpImage,1),1)';
allImage = nan([size(tmpImage,1)*size(tmpImage,2),length(x)]);
progBar = ProgressBar(length(x),'calculating distance...');
for i = 1:length(x)
    % get matrix elements, based on sigma
    %distMat(:,:,i) = sqrt( (yMat - tmpY(i)).^2 + (xMat - tmpX(i)).^2);
    distMat = sqrt( (yMat - tmpY(i)).^2 + (xMat - tmpX(i)).^2);
    allImage(distMat< sig(i)*2,i) = in(i);
    
    
    %%%% CHANGE sig(i) TO A SCALAR VALUE, as voxels already are scaled by eccentricity %%%%%   
    
    progBar(i);
end
newImage = reshape(allImage,[size(tmpImage,1),size(tmpImage,2),length(x)]);
meanImage = nanmean(newImage,3);
H = fspecial('gaussian',[2 2],0.5);
blurMat = imfilter(meanImage,H);
blurMat(blurMat == 1) = nan;
pcolor(blurMat);
shading flat;
%imagesc(blurMat,[0 0.5]);
%colormap('hot');
hold on;
% Plot circles
for i = 1:length(cLines)
    plot(cX(i,:),cY(i,:),'k:','LineWidth',1);
    if cLines(i) ~= 0
        [tX,tY] = pol2cart(deg2rad(15),cLines(i)*pscale);
        text(-tX + max(cX(:))/2,tY + max(cY(:))/2,num2str(round(10^cLines(i))), 'FontSize', 20,...
            'HorizontalAlignment','center','VerticalAlignment','middle');
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
%% Set lims
% logaxLim = log10(axLim);
% xlim([-logaxLim logaxLim]);
% ylim([-logaxLim logaxLim]);
%% Plot data
switch plotType
    case 'scatterPlot'
        in(in>caxLims(2)) = caxLims(2);
        scatter(x,y,sigpix,in,'filled');
    case 'fillPlot'
        
        
        H = fspecial('gaussian',[10 10],2);
blurMat = imfilter(tmpImage,H);
figure;
subplot(1,2,1);
imshow(tmpImage);
subplot(1,2,2);
imshow(blurMat);
        
        
        for i = 1:numel(gridX)
            % find closest points
            
        end
end

%% Plot legends
axes('position',[0.85 0.1 0.1 0.8]) ; % inset
pY = 10*(1:5);
pX = ones(size(pY));
plotsigs = linspace(ceil(min(sig)),ceil(max(sig)),5);
plotpix = linspace(min(sigpix),max(sigpix),5);
scatter(pX,pY,plotpix,'filled');
%th = title('Sigma size','FontSize',20);
%P = get(th,'Position');
%set(th,'Position',[P(1),P(2),P(3)]);
tX = 0*pX;
tY = pY;
for i = 1:length(tX)
    text(tX(i),tY(i),sprintf('%0.2f',plotsigs(i)),'FontSize',20,...
        'HorizontalAlignment','center','VerticalAlignment','middle');
end
axis off;
% Colorbar
h=colorbar;
set(h,'fontsize',20);
%yl = ylabel(h,'Data value','FontSize',20,'rot',-90);
%yP = get(yl,'Position');
%set(yl,'Position',[yP(1),yP(2),yP(3)]);
caxis(caxLims)