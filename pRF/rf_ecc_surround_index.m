function [I] = rf_ecc_surround_index(E,visual_area)

% Calculates the receptive field suppressive index for a given eccentricity
%   (E), for a visual area (V1, V2, V3).
%
%   Usage:
%   [D] = rf_ecc_surround_index(E,visual_area)
%
%       I = index
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
    P = [-0.06189 0.85946];
elseif strcmp(visual_area,'V2')
    P = [-0.02556 0.65911];
elseif strcmp(visual_area,'V3')
    P = [0.04855 0.25753];
end
I = polyval(P,E);
