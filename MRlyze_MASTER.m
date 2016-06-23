%% Template script for analyzing MRI data using MRlyze
%
%   Required software:
%       Freesurfer, FSL
%       IMPORTANT:
%           start matlab from terminal, so environmental variables, 
%               libraries are set correctly
%           add $FREESURFER_HOME/matlab to your matlab path
%
%   Directory structure:
%       data_directory -> project directory -> subject directory ->
%           session directory -> dicom directory
%
%           e.g. ~/data/Retinotopy/ASB/10012014/DICOMS
%
%       If physiological measures are collected, using the following 
%           directory structure:
%
%       data_directory -> project directory -> subject directory ->
%           session directory -> physio directory
%
%           e.g. ~/data/Retinotopy/ASB/10012014/PulseOx
%
%   I recommend creating a project specific master file in a separate 
%   directory. For example, copy this MASTER.m file to a project 
%   specific folder:
%
%   /User/Shared/Matlab/<project_name>/<project_name>_MASTER.m
%
%   Written by Andrew S Bock Jun 2016

