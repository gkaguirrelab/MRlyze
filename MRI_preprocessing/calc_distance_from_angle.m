function [d] = calc_distance_from_angle(a,f)

% Calculates the distance (d) given an angle (a) and focal length (f)
%
%   Usage:
%   [d] = calc_distance_from_angle(a,f)
%
%   note: assumes angle is in radians
%
%   Written by Andrew S Bock Sep 2015
%%
d = tan(a/2) * (2*f);