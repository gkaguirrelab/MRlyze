function [D] = rf_ecc(E,visual_area)

% Calculates the receptive field size in degrees of visual angle for a
%   given eccentricity (E), for a visual area (V1, V2, V3).
%
%   Usage:
%   [D] = rf_ecc(E,visual_area)
%
%       D = degrees of visual angle
%       E = eccentricty in degrees
%       visual_area = 'V1' <default> 'V2' or 'V3'
%
%   Equations were created from data in figure 6 in:
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
%% Calculate sigma at given eccentricty
if strcmp(visual_area,'V1')
    P = [0.43451 0.83335];
elseif strcmp(visual_area,'V2')
    P = [0.59913 0.91096];
elseif strcmp(visual_area,'V3')
    P = [0.85831 1.54683];
end
D = polyval(P,E); % values in FWHM
D = D./(2*sqrt(2*log(2))); % convert to sigma
%% Below are from Dumoulin and Wandell (2008)
% if strcmp(visual_area,'V1')
%     P = [0.05770 0.32159];
% elseif strcmp(visual_area,'V2')
%     P = [0.09776 0.44012];
% elseif strcmp(visual_area,'V3')
%     P = [0.15096 0.99046];
% end
%D = polyval(P,E);

