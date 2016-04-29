function [commandc] = feat_HCLV_Photo_7T(outFile,funcVol,anatVol,EVs,WA)

% Write a .fsf file for first-level feat for the HCLV_Photo_7T protocol
%
%   Usage:
%   feat_HCLV_Photo_7T(outFile,funcVol,anatVol,EVs,WA)
%
%   Inputs:
%   outFile     = name of output .fsf file (fullpath)
%   funcVol     = name of input functional volume (fullpath)
%   anatVol     = name of 'standard' image for registration (fullpath)
%   EVs         = structure containing EV text files for each condition (fullpath)
%   WA          = 1 - wrap around block at beginning of run; 0 - no wrap around
%
%   Written by Andrew S Bock Apr 2016

%% Set defaults
tmpDir = which('feat_HCLV_Photo_7T');
design_dir = fileparts(tmpDir);
%% Load in template
if WA
    design_file = fullfile(design_dir,'first-level-template_HCLV_Photo_7T_WA.fsf');
else
    design_file = fullfile(design_dir,'first-level-template_HCLV_Photo_7T.fsf');
end
%% Load functional volume
tmp = load_nifti(funcVol);
%% Set design values
DESIGN.TR = num2str(tmp.pixdim(5)/1000); % TR is in msec, convert to sec
if tmp.pixdim(5) < 100 % use 100, in case very short TR is used (i.e. multi-band)
    error('TR is not in msec');
end
DESIGN.VOLS = num2str(tmp.dim(5));
DESIGN.STANDARD = anatVol;
DESIGN.TOTAL_VOXELS = num2str(tmp.dim(2)*tmp.dim(3)*tmp.dim(4)*tmp.dim(5));
DESIGN.FEAT_DIR = funcVol;
for i = 1:length(EVs)
    eval(['DESIGN.EV' num2str(i) ' = ' EVs{i}]);
end
fin = fopen(design_file,'rt');
fout = fopen(outFile,'wt');
fields = fieldnames(DESIGN);
while(~feof(fin))
    s = fgetl(fin);
    for f = 1:length(fields)
        s = strrep(s,['DESIGN_' fields{f}],DESIGN.(fields{f}));
    end
    fprintf(fout,'%s\n',s);
    %disp(s)
end
fclose(fin);
% Diplay the feat command
commandc = ['feat ' outFile];
disp(commandc);