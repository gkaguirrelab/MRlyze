function [outImage] = plotVisual(params)

% Plots values from an input volume in visual field coordinates
%
%   Usage:
%   [outImage] = plotVisual(params)
%
%   See also:
%   convert_image2surf
%
%   Written by Andrew S Bock Mar 2016
%   Updates Dec 2016 based on Bosco Tjan's `vfmap2` function

%% set defaults
if ~isfield(params,'maxextent')
    params.maxextent        = 15;
end
if ~isfield(params,'visualArea')
    params.visualArea        = 'V1';
end
if ~isfield(params,'gridsize')
    params.gridsize          = 501; % should be an odd number to include zero
end
if ~isfield(params,'intpmethod')
    params.intpmethod       = 'natural'; % interpolation method for griddata
end
if ~isfield(params,'ringstep')
    params.ringstep         = 5;     % deg
end
if ~isfield(params,'textgap')
    params.textgap          = 0.15;	% proportion of ringstep
end
if ~isfield(params,'textsize')
    params.textsize         = 12;
end
if ~isfield(params,'textcolor')
    params.textcolor        = [1,1,1];
end
if ~isfield(params,'bgcolor')
    params.bgcolor          = [0.2 0.2 0.2];
end
if ~isfield(params,'caxis')
    params.caxis            = [-8 8]; % color range
end
%% Make custom color map
mymap                       = zeros([64 3]);
mymap(1:32,1)               = linspace(1,0,32);
mymap(33:end,2)             = linspace(0,1,32);
mymap(1:8,2)                = linspace(0.5,0,8);
mymap(1:8,3)                = linspace(0.5,0,8);
mymap(end-7:end,3)          = linspace(0,0.5,8);
mymap(end-7:end,1)          = linspace(0,0.5,8);
params.cmap                 = mymap;
%% Load in volumes
% Eccentricity 
inData = load_nifti(params.lhEcc);
ecc = inData.vol;
inData = load_nifti(params.rhEcc);
ecc = [ecc;inData.vol];
ecc(isnan(ecc)) = 0;
% Polar angle 
inData = load_nifti(params.lhPol);
pol = inData.vol;
inData = load_nifti(params.rhPol);
pol = [pol;inData.vol];
pol(isnan(pol)) = 0;
% Areas
inData = load_nifti(params.lhAreas);
areas = inData.vol;
inData = load_nifti(params.rhAreas);
areas = [areas;inData.vol];
areas(isnan(areas)) = 0;
% Values
inData = load_nifti(params.lhVals);
vals = inData.vol;
inData = load_nifti(params.rhVals);
vals = [vals;inData.vol];
vals(isnan(vals)) = 0;
%% Threshold based on visual areas
switch params.visualArea
    case 'V1'
        goodInd = abs(areas)==1;
    case 'V2'
        goodInd = abs(areas)==2;
    case 'V3'
        goodInd = abs(areas)==3;
end
ecc(~goodInd) = [];
pol(~goodInd) = [];
areas(~goodInd) = [];
vals(~goodInd) = [];
%% Convert polar to cartesian
[x,y] = pol2cart(pol,ecc);
% flip y
y = -y;
%% Plot the data points in visual space
figure('units','normalized','position',[0 0 1 1]);
polarplot(pol,ecc,'ko','MarkerFaceColor','k','MarkerSize',1);
rlim([0 15]);
ax = gca;
ax.Color = [0.5 0.5 0.5];
ax.RColor = [1 1 1];
ax.GridColor = [1 1 1];
ax.GridAlpha = 1;
%% Plot the values in visual space
figure('units','normalized','position',[0 0 1 1]);
onhold = ishold;
% set the axis
caxis(params.caxis);
xlim([-params.maxextent params.maxextent]*1.1); % add 1% to accommodate thick grid lines
ylim([-params.maxextent params.maxextent]*1.1);
axis equal
hold on

[x,y] = pol2cart(pol,ecc);
y = -y;
[X,Y] = meshgrid(linspace(-params.maxextent,params.maxextent,params.gridsize));
vq = griddata(x,y,vals,X,Y,params.intpmethod);
%vq(abs(vq)<params.toosmall) = 0;
vq(X.^2+Y.^2>params.maxextent^2) = nan;
% pcolor(X,Y,vq)
surf(X,Y,zeros(size(X)),vq)
zlim([min(vq(:)) max(vq(:))])
shading flat
colormap(params.cmap);
grid off
axis off
axis image
% plot the grid
lim = floor(params.maxextent/params.ringstep);
gap = params.ringstep*params.textgap;
text(gap,gap, '0','FontSize',params.textsize,'Color',params.textcolor);
for i=1:lim
    h=mycircle([0 0],i*params.ringstep); h.EdgeColor = params.textcolor;
    text(0+gap,i*params.ringstep+gap, num2str(i*params.ringstep),'FontSize',params.textsize,'Color',params.textcolor);
end
line([-lim*params.ringstep lim*params.ringstep], [0 0], 'color', params.textcolor);
line([0 0], [-lim*params.ringstep lim*params.ringstep], 'color', params.textcolor);
h=colorbar; h.Color = params.textcolor;

if ~onhold
    hold off
end
H = gcf;
H.Color = params.bgcolor;