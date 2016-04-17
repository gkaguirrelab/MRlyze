function [varexp,params,dists,newx,newy,newz] = plot_template_fits(tdir,template,hemi,fitType,makeplot,crange)

% Makes a 3D scatter plot of pRF template fits, following
% 'regress_template'
%
%   Usage:
%   [varexp,params,dists,newx,newy,newz] = plot_template_fits(session_dir,hemi,makeplot,crange)
%
%   Inputs:
%   session_dir = session directory
%   template = 'coarse' or 'fine'
%   hemi = 'lh' or 'rh'
%   makeplot = 0 <default> or 1 % makes a 3D scatter plot
%   crange = color range for colorbar <default> [0 0.05]
%
%   Written by Andrew S Bock Sep 2015

%% Set defaults
if ~exist('fitType','var')
    fitType = 'V1';
end
if ~exist('makeplot','var')
    makeplot = 0;
end
if ~exist('crange','var')
    crange = [0 0.05];
end
%% Get the template parameters and associated variance explained.
[varexp,params] = find_best_template(template,tdir,hemi,[],[],[],fitType);
for i = 1:length(varexp)
    x(i) = params(i).FCx;
    y(i) = params(i).FCy;
    z(i) = params(i).psi;
end
%% Determine the 'distance' from the best template (highest variance explained)
for i = 1:length(varexp)
    newx = x - x(1);
    newy = y - y(1);
    newz = z - z(1);
end
dists = sqrt( (x - x(1)).^2 + (y - y(1)).^2 + (z - z(1)).^2 );
%% Plot 3D scatter plot
if makeplot
    dotSize = 1000;
    figure;scatter3(newx,newy,newz,dotSize,varexp,'filled');
    colorbar;caxis(crange);
end