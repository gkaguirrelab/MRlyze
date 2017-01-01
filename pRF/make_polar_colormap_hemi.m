function [mycolormap] = make_polar_colormap_hemi(mapres,show)

%   Creates a matrxi useful for plotting pRF and ccRF polar angle maps
%
%   Usage:
%   [mycolormap] = make_polar_colormap(show)
%
%   defaults:
%   show = 0; do not plot the resulting polar angle colormap
%
%   Written by Andrew S Bock Oct 2014

%% Set up defaults
if ~exist('show','var')
    show = 0;
end
%% Create colormap
%mycolormap = blue_green_red;
mycolormap = .8*ones(mapres(3),3);
% Red to Yellow (1,0,0) -> (1,0.85,0)
mycolormap(2*mapres(3)/8+1:3*mapres(3)/8,1) = ones(mapres(3)/8,1);
mycolormap(2*mapres(3)/8+1:3*mapres(3)/8,2) = linspace(0,0.85,mapres(3)/8);
mycolormap(2*mapres(3)/8+1:3*mapres(3)/8,3) = zeros(mapres(3)/8,1);
% Yellow to Green (1,0.85,0) -> (0,0.75,0)
mycolormap(3*mapres(3)/8+1:4*mapres(3)/8,1) = linspace(1,0,mapres(3)/8);
mycolormap(3*mapres(3)/8+1:4*mapres(3)/8,2) = linspace(0.85,0.75,mapres(3)/8);
mycolormap(3*mapres(3)/8+1:4*mapres(3)/8,3) = zeros(mapres(3)/8,1);
% Green to Cyan (0,0.75,0) -> (0,0.85,1)
mycolormap(4*mapres(3)/8+1:5*mapres(3)/8,1) = zeros(mapres(3)/8,1);
mycolormap(4*mapres(3)/8+1:5*mapres(3)/8,2) = linspace(0.75,0.85,mapres(3)/8);
mycolormap(4*mapres(3)/8+1:5*mapres(3)/8,3) = linspace(0,1,mapres(3)/8);
% Cyan to Blue (0,0.85,1) -> (0,0,1)
mycolormap(5*mapres(3)/8+1:6*mapres(3)/8,1) = zeros(mapres(3)/8,1);
mycolormap(5*mapres(3)/8+1:6*mapres(3)/8,2) = linspace(0.85,0,mapres(3)/8);
mycolormap(5*mapres(3)/8+1:6*mapres(3)/8,3) = ones(mapres(3)/8,1);
% Flip so visual field is reversed
mycolormap = flipud(mycolormap);
%% Plot resulting colormap
if show
    figure;
    x = linspace(-1,1,mapres(3));
    [xx,yy] = meshgrid(x);
    [th,r] = cart2pol(xx,yy);
    img = th;
    img = img/nanmax(img(:))*mapres(2);
    img(r>1) = -3; % Set values outside of circle to black
    imagesc(x,x,img);
    colorbar
    axis tight equal off
    set(gca,'ydir','normal')
    colormap(mycolormap)
end