function [a] = calc_visual_angle(d,f)

% Calculates the visual angle(s) (a) given a distance(s) (d) and focal 
%   length (f). Output (a) is in radians.
%
%   Usage:
%   [a] = calc_visual_angle(d,f)
%
%   Written by Andrew S Bock Sep 2015

%%
a = 2.*atan(d./(2.*f));