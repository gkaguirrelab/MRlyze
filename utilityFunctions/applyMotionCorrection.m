function applyMotionCorrection(inFile,outFile,motDir,refvol)

% Applys motion correction using MAT files (FSL style)
%
%   Usage:
%   applyMotionCorrection(inFile,outFile,motDir,refNum)
%
%   This is typically used to apply motion parameters (filtered at the 
%   Nyquist frequency), created using 'filterMotionMatFiles'
%
%   Note:
%   refNum = number of reference TR (starts at 1, not 0!)
%
%   Based on scripts by Marta Vidorreta Díaz de Cerio
%   Written by Andrew S Bock Apr 2016

%% set defaults
if ~exist('refvol','var')
   refvol = 1; % defaults is first TR 
end
%% Apply motion correction
system(['fslsplit ' inFile ' grot']);
listvols = '';
nTRs = listdir('./grot*.nii.gz','files');
n_scans = length(nTRs);
for iscan = 1:n_scans
    commandc = sprintf(['applywarp -i grot%04d -o grot%04d --premat=' ...
        fullfile(motDir,'MAT_%04d') ' -r grot%04d --rel'],iscan-1,iscan-1,iscan,refvol-1);
    disp(commandc);
    system(commandc);
    listvols = [listvols ' ' sprintf('grot%04d', iscan-1)];
end
commandc = ['fslmerge -t ' outFile ' ' listvols];
disp(commandc);
system(commandc);
commandc = 'rm -f grot*';
system(commandc);