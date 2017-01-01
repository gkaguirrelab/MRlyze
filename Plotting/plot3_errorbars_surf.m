function [h]=plot3_errorbars_surf(x,y,z,e,aLims,aLabels,MSize,LWidth,GSize,FSize)
% Plots 3d data using the plot3 and surf functions. Adds vertical errorbars 
%   to each point
%
%   Usage:
%   [h]=plot3_errorbars_surf(x,y,z,e,aLims,aLabels,MSize,LWidth,GSize)
%
%   Inputs:
%   x = xaxis
%   y = yaxis
%   z = zaxis
%   e = error value
%   aLims = axes limits <default = {[-15 15],[-15 15],[0 0.05]}>
%   aLabels = axes labels <default = {'x-axis' 'y-axis' 'z-axis'}>
%   MSize = error marker size <default = 10>
%   LWidth = error line width <default = 1>
%   GSize = grid size for surface mesh <default = 0.1>
%
%   Taken from:
%   http://code.izzid.com/2007/08/19/How-to-make-a-3D-plot-with-errorbars-in-matlab.html
%   Article created: Aug 19, 2007
%   Article by: Jeremiah Faith
%
%   Updated by Andrew S Bock Sep 2015

%% set defaults
if ~exist('aLims','var')
    aLims = {[-15 15],[-15 15],[0 0.05]};
end
if ~exist('aLabels','var')
    aLabels = {'x-axis' 'y-axis' 'z-axis'};
end
if ~exist('MSize','var')
    MSize = 10;
end
if ~exist('LWidth','var')
    LWidth = 1;
end
if ~exist('GSize','var')
    GSize = 0.1; % Grid size
end
if ~exist('FSize','var')
    FSize = 20; % Font Size
end
%% Make the plot
% create the standard 3d scatterplot
hold off;
ph=plot3(x, y, z, '.k','LineWidth',0.01,'MarkerSize', MSize);
% looks better with large points
hold on
% now draw the vertical errorbar for each point
for i=1:length(x)
    xV = [x(i); x(i)];
    yV = [y(i); y(i)];
    zMin = z(i) + e(i);
    zMax = z(i) - e(i);
    zV = [zMin, zMax];
    % draw vertical error bar
    eh=plot3(xV, yV, zV, '-k');
    set(eh, 'LineWidth', LWidth,'MarkerSize', MSize);
end
% now we want to fit a surface to our data
% the  0.25 and 0.1 define the density of the fit surface
% adjust them to your liking
tt1=[floor(min(min(x))):GSize:max(max(x))];
tt2=[floor(min(min(y))):GSize:max(max(y))];
% prepare for fitting the surface
[xg,yg]=meshgrid(tt1,tt2);
% fit the surface to the data;
% matlab has several choices for the fit;  below is "linear"
zg=griddata(x,y,z,xg,yg,'linear');
% draw the surf on our plot
sh = surf(xg,yg,zg,'EdgeColor','none','LineStyle','none');
xlabel(aLabels{1},'FontSize',FSize);
ylabel(aLabels{2},'FontSize',FSize);
zlabel(aLabels{3},'FontSize',FSize);
% deal with axes
xlim(aLims{1});
ylim(aLims{2});
zlim(aLims{3});
axis square;
grid on;