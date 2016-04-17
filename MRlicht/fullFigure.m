function [h,screensize] = fullFigure
% Creates a full screen figure window
%
%   Usage:
%   h = fullFigure;
%
%   Written by Andrew S Bock Mar 2016

%% Create figure
screensize = get(0,'Screensize');
h = figure('units','normalized','position',[0 0 1 1]);