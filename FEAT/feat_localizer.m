function feat_localizer(session_dir,runNum,func,run_feat,design_file)

% Creates a .fsf file in the 'runNum' bold directory for a given
% session directory. For use in localization of areas.
%
%   Usage:
%   feat_localizer(session_dir,runNum,func,run_feat,design_file)
%
%   Defaults:
%     func = 's2.rf.tf'; % default functional input
%     run_feat = 0; % default will NOT run feat, just create .fsf file
%     design_file = fullfile(design_dir,'feat_localizer_template.fsf');
%
%   Written by Andrew S Bock Nov 2015

%% Set default parameters
if ~exist('func','var')
    func = 's2.rf.tf';
end
if ~exist('run_feat','var')
    run_feat = 0; % don't run FEAT, e.g. if on a cluster
end
if ~exist('design_file','var')
    design_dir = which('feat_TTF');
    design_dir = fileparts(design_dir);
    design_file = fullfile(design_dir,'feat_localizer_template.fsf');
end
%% Add to log
SaveLogInfo(session_dir,mfilename,session_dir,runNum,func,run_feat,design_file);
%% Find bold run directories
d = find_bold(session_dir);

%% Find structural file
if exist(fullfile(session_dir,'MP2RAGE'),'dir')
    out_file = fullfile(session_dir,'MP2RAGE','004','MP2RAGE_brain.nii.gz'); % brain extracted output file
    %skull_strip(session_dir,subject_name);
    DESIGN.STRUCT = out_file;
elseif exist(fullfile(session_dir,'MPRAGE'),'dir')
    out_file = fullfile(session_dir,'MPRAGE','001','MPRAGE_brain.nii.gz'); % brain extracted output file
    %skull_strip(session_dir, subject_name);
    DESIGN.STRUCT = out_file;
end
%% Run FEAT
tmp = load_nifti(fullfile(session_dir,d{runNum},[func '.nii.gz']));
DESIGN.TR = num2str(tmp.pixdim(5)/1000); % TR is in msec, convert to sec
if tmp.pixdim(5) < 100 % use 100, in case very short TR is used (i.e. multi-band)
    error('TR is not in msec');
end
DESIGN.VOLS = num2str(tmp.dim(5));
DESIGN.STANDARD = DESIGN.STRUCT;
DESIGN.TOTAL_VOXELS = num2str(tmp.dim(2)*tmp.dim(3)*tmp.dim(4)*tmp.dim(5));
DESIGN.FEAT_DIR = fullfile(session_dir,d{runNum},[func '.nii.gz']);
DESIGN.CONDITION_1 = fullfile(session_dir,d{runNum},'condition_1.txt');
DESIGN.CONDITION_2 = fullfile(session_dir,d{runNum},'condition_2.txt');
fin = fopen(design_file,'rt');
outname = [func '.feat_localizer.fsf'];
fout = fopen(fullfile(session_dir,d{runNum},outname),'wt');
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
disp('Running:')
commandc = ['feat ' fullfile(session_dir,d{runNum},outname)];
disp(commandc);
% Run feat
if run_feat
    [~,~] = system([commandc ' &']);
end
% Also save this out in a file.
fid = fopen(fullfile(session_dir, 'feat_localizer_scripts'), 'a');
fprintf(fid, '%s\n',commandc);
fclose(fid);