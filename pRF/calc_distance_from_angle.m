function d = calc_distance_from_angle(a,f)

% Calculates the distance(s) (d) given angle(s) (a) and focal length(s) (f). 
%   Input angle(s) in radians.
%
%   Usage:
%   [d] = calc_distance_from_angle(a,f)
%
%   Written by Andrew S Bock Sep 2015

%%
d = tan(a./2) .* (2.*f);