function [M] = cortical_mag(E,visual_area)

% Calculates the mm/deg cortical magnification for a given eccentricity (E),
%   for a visual area (V1, V2, or V3). 
%
%   Usage:
%   [M] = cortical_mag(E,visual_area)
%
%       M = mm/deg
%       E = eccentricty in degrees
%       visual_area = 'V1' <default> 'V2' or 'V3'
%
%   Values are calculated using equations from:
%
%   Dougherty, R. F., Koch, V. M., Brewer, A. A., Fischer, B., Modersitzki,
%       J., & Wandell, B. A. (2003). Visual field representations and
%       locations of visual areas V1/2/3 in human visual cortex. Journal of
%       Vision, 3(10), 1.
%
%   Note (from that paper): 
%       Because we did not measure eccentricities closer than 2°, the e2 
%       estimates are not robust. Also, the variance across individuals is 
%       large, especially at very central locations (see Figure 4). 
%
%   Written by Andrew S Bock Jul 2015

%% Set defaults
if ~exist('visual_area','var')
    visual_area = 'V1';
end
%% Calculate mm/deg at given eccentricty
if strcmp(visual_area,'V1')
    A = 29.2; % millimeters
    e2 = 3.67; % degrees
elseif strcmp(visual_area,'V2')
    A = 22.8; % millimeters
    e2 = 2.54; % degrees
elseif strcmp(visual_area,'V3')
    A = 19.4; % millimeters
    e2 = 2.69; % degrees
end
M = A ./ (E + e2);
