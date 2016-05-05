function FIR_first_level_feat(outFile,funcVol,anatVol,EVs,condition)

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
matlab_path = '/data/jag/MELA/Matlab';
design_dir = fullfile(matlab_path,'MelanopsinMR/FEAT_templates');
%% Load in template
switch condition
    case 'MelPulses_400pct'
        design_file = fullfile(design_dir,'FIR_MaxMel_template.fsf');
    case 'LMSPulses_400pct'
        design_file = fullfile(design_dir,'FIR_MaxLMS_template.fsf');
    case 'SplatterControl'
        design_file = fullfile(design_dir,'FIR_SplatterControl_template.fsf');
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
    eval(['DESIGN.EV' num2str(i,'%02d') ' = ''' EVs{i} ''';']);
end
disp(DESIGN);
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
