function [mycolormap] = make_colormap(color_components)

%   Creates a colormap with custom transition colors
%
%   Usage: 
%   [mycolormap] = make_polar_colormap(color_components,show)
%   color_components: k x 3 matrix with each color component in one row
%
%   defaults:
%   show = 0; do not plot the resulting polar angle colormap
%
%   Written by Andrew S Bock Oct 2014

%% Set up defaults
if ~exist('color_components','var')
    color_components = [1 0 0;0 0 1];
end
%% Create colormap

% Check if dimensions are inverted and correct if so
if and(size(color_components,1)==3,size(color_components,2)~=3)
    color_components = color_components';
end
% Check if color scale is 0-255 and rescale to 0-1 if so
if max(color_components(:))>1
    color_components = color_components/255;
end
numColors = size(color_components,1);
mycolormap = zeros(100*(numColors-1),3);

% Create color map
for c=1:(numColors-1)
    mycolormap(((c-1)*100+1):(c*100),:) = repmat(color_components(c,:),100,1) + (0:0.01:0.99)'*(color_components(c+1,:)-color_components(c,:));
end
% Add final color to the last row
mycolormap(end+1,:) = color_components(end,:);