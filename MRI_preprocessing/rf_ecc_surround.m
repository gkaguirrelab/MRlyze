function [D] = rf_ecc_surround(E,visual_area)

% Calculates the receptive field surround size in degrees of visual angle 
%   for a given eccentricity (E), for a visual area (V1, V2, V3).
%
%   Usage:
%   [D] = rf_ecc_surround(E,visual_area)
%
%       D = degrees of visual angle
%       E = eccentricty in degrees
%       visual_area = 'V1' <default> 'V2' or 'V3'
%
%   Equations were created from data in figure 7 in:
%
%   Zuiderbaan, W., Harvey, B. M., & Dumoulin, S. O. (2012). Modeling 
%       center?surround configurations in population receptive fields using 
%       fMRI. Journal of vision, 12(3), 10.
%
%   Average fMRI-pRF estimate values were extracted using WebPlotDigitizer
%   (http://www.arohatgi.info/WebPlotDigitizer/)
%
%   Written by Andrew S Bock Jul 2015

%% Set defaults
if ~exist('visual_area','var')
    visual_area = 'V1';
end
%% Calculate mm/deg at given eccentricty
if strcmp(visual_area,'V1')
    P = [0.88403 7.85043];
elseif strcmp(visual_area,'V2')
    P = [1.30916 7.78662];
elseif strcmp(visual_area,'V3')
    P = [1.76160 10.06623];
end
D = polyval(P,E); % values in FWHM
D = D./(2*sqrt(2*log(2))); % convert to sigma
