function mri_surf2surf(srcsubject,trgsubject,sval,tval,hemi)

% Uses Freesurfer's 'mri_surf2surf' command to project from a surface
% volume from one surface to another
%
%   Usage:
%   mri_surf2surf(srcsubject,trgsubject,sval,tval,hemi)
%
%   Written by Andrew S Bock Nov 2015

%% Project surface
system(['mri_surf2surf --srcsubject ' srcsubject ' --trgsubject ' trgsubject ...
    ' --sval ' sval ' --tval ' tval ' --hemi ' hemi]);
