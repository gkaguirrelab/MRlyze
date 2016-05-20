function plotVisual(inVol,eccVol,polVol,sigVol,roiInd,roiName,axLim)

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
circPol = 0:pi/50:2*pi; % Polar angle sampling
circEcc = ones(size(circPol)); % Eccentricity sampling
cLines = linspace(0,log10(axLim),6); % Eccentricity lines
rLines = linspace(0,(2*pi) - (2*pi)/12,12); % Polar angle lines (spokes)
pLabels = cLines(end) + (cLines(end) - cLines(end-1))/2; % Polar angle label eccentricity
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
% Convert sig (deg) to sig (pixels)
sigpix = angle2pix(display,sig);
% load input volume
in = load_nifti(inVol);
in = in.vol(roiInd);
%% Convert polar to cartesian
[x,y] = pol2cart(pol,logecc);
% flip y
y = -y;
%% Plot circles and spokes
% Create circles
for i = 1:length(cLines)
    [cX(i,:),cY(i,:)] = pol2cart(circPol,cLines(i)*circEcc);
end
% Grey background
patch('XData',cX(end,:),'YData',cY(end,:),...
    'FaceColor',[1 1 1],'HandleVisibility', 'off');
hold on;
% Plot circles
for i = 1:length(cLines)
    plot(cX(i,:),cY(i,:),'k:','LineWidth',1);
end
% Plot spokes
for i = 1:length(rLines)
    [X,Y] = pol2cart(rLines(i),cLines(end));
    plot([X -X],[Y -Y],'k:','LineWidth',1);
    [tX,tY] = pol2cart(rLines(i),pLabels);
    text(tX,tY,num2str(rad2deg( rLines(i) )) , 'FontSize', 20,...
        'HorizontalAlignment','center','VerticalAlignment','middle');
    text(-tX,-tY,num2str(rad2deg( mod(pi+rLines(i),2*pi) )) ,'FontSize', 20,...
        'HorizontalAlignment','center','VerticalAlignment','middle');
end
% Set lims
logaxLim = log10(axLim);
xlim([-logaxLim logaxLim]);
ylim([-logaxLim logaxLim]);
axis square;
%% Plot data
scatter(x,y,sigpix,in,'filled');

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
caxis([min(in) max(in)])