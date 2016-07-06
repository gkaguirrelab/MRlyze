function FIR_first_level_feat(outFile,funcVol,anatVol,EVs,condition)
% FIR_first_level_feat(outFile,funcVol,anatVol,EVs,condition)
%
% Writes a .fsf file for first-level feat for the MelanopsinMR Protocols
% based of pregenerated template files.
%
%
% Input arguments:
% ================
%
%   outFile : name of the fsf output file;
%   funcVol : functional volume to use in the analysis;
%   anatVol : anatomical volume to use in the analysis;
%   EVs : path to the EV files to use in the analysis.
%   condition : name of the condition for the analysis. The fsf template
%   wil be selected accordingly.
% 
% Usage:
% ======
%   This function is better used in a script , looping through different
%   bold runs. An example is shown in
%   'MelanopsinMR/analysis/MelanopsinMR_FeatStatAnalysis.m' lines 125-179.
% 
%
%
% 7/2/2016  gf, ms      Written and commented.

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
    case 'MaxMelCRF'
        design_file = fullfile(design_dir,'MaxMelCRF_template.fsf');
    case 'MaxLMSCRF'
        design_file = fullfile(design_dir,'MaxLMSCRF_template.fsf');
    otherwise
        error('No template for this condition. Try using: ''MelPulses_400pct'', ''LMSPulses_400pct'', ''SplatterControl'', ''MaxMelCRF'' or ''MaxLMSCRF''');
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
    eval(['DESIGN.EV' num2str(i,'%03d') ' = ''' EVs{i} ''';']);
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
