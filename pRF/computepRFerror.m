function [degerror] = computepRFerror(inEcc,inPol,tempEcc,tempPol,verts)

% Computes the error in degrees visual angle between input eccentricity and
% polar angle data and a template
%
%   Usage:
%   [error] = computepRFerror(inEcc,inPol,tempEcc,tempPol,verts)
%
%   Inputs:
%   inEcc       = input surface file of eccentricity values
%   inPol       = input surface file of polar angle values
%   tempEcc     = template surface file of eccentricity values
%   tempPol     = template surface file of polar angle values
%   verts       = logical vector specifying the vertices of interest
%
%   Output:
%   degerror    = difference (deg) between template and input for each vertex 
%
%   Written by Andrew S Bock Jun 2016

%% Load in data
eccIn       = load_nifti(inEcc);
polIn       = load_nifti(inPol);
eccTemp     = load_nifti(tempEcc);
polTemp     = load_nifti(tempPol);
%% Computer error in visual angle
[inX,inY] = pol2cart(polIn.vol(verts),eccIn.vol(verts));
[tempX,tempY] = pol2cart(polTemp.vol(verts),eccTemp.vol(verts));
degerror = sqrt((inX - tempX).^2 + (inY - tempY).^2);