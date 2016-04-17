function [outMat] = convertMotion2Mat(motParams)

% Takes in six motion params (3 rotations, 3 translations) and converts to
% a transformation matrix (FSL style).
%
%   Usage:
%   [outMat] = convertMotion2Mat(motParams)
%
%   Based on scripts by Marta Vidorreta Díaz de Cerio
%   Written by Andrew S Bock Apr 2016

%% Pull out individual parameters
xalpha = motParams(1);
yalpha = motParams(2);
zalpha = motParams(3);
tx     = motParams(4) ;
ty     = motParams(5) ;
tz     = motParams(6) ;
%% Rotation X
Rx = [...
    1               0               0               0 ;
    0               cos(xalpha)     sin(xalpha)     0 ;
    0               -sin(xalpha)    cos(xalpha)     0 ;
    0               0               0               1];
%% Rotation Y
Ry = [...
    cos(yalpha)     0               -sin(yalpha)    0 ;
    0               1               0               0 ;
    sin(yalpha)     0               cos(yalpha)     0 ;
    0               0               0               1];
%% Rotation Z
Rz = [...
    cos(zalpha)     sin(zalpha)     0               0 ;
    -sin(zalpha)    cos(zalpha)     0               0 ;
    0               0               1               0 ;
    0               0               0               1];
%% Translations
T = [...
    0               0               0               tx;
    0               0               0               ty;
    0               0               0               tz;
    0               0               0               0];
%% Create output matrix
M = Rx * Ry * Rz + T;
outMat = M';