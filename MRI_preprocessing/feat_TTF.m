function feat_TTF(session_dir,func,run_feat,design_file)

% Creates a 'feat_TTF.fsf' file in each bold directory for a given
% session directory. For use in subsequent create of TTFs.
%
%   Usage:
%   feat_TTF(session_dir,func,run_feat,design_file)
%
%   Defaults:
%     func = 'dbrf.tf'; % default functional input
%     run_feat = 0; % default will NOT run feat
%     design_file = 'feat_TTF_template.fsf';
%
%   Written by Andrew S Bock Oct 2015

%% Set default parameters
if ~exist('func','var')
    func = 'dbrf.tf';
end
if ~exist('run_feat','var')
    run_feat = 0; % run FEAT from matlab.  Alternative is to run from the terminal, e.g. if on a cluster
end
if ~exist('design_file','var')
    design_dir = which('feat_TTF');
    design_dir = fileparts(design_dir);
    design_file = fullfile(design_dir,'feat_TTF_template.fsf');
end
%% Add to log
SaveLogInfo(session_dir,mfilename,session_dir,func,run_feat,design_file);
%% Find bold run directories
d = find_bold(session_dir);
nruns = length(d);
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
for r = 1:nruns
    tmp = load_nifti(fullfile(session_dir,d{r},[func '.nii.gz']));
    DESIGN.TR = num2str(tmp.pixdim(5)/1000); % TR is in msec, convert to sec
    if tmp.pixdim(5) < 100 % use 100, in case very short TR is used (i.e. multi-band)
        error('TR is not in msec');
    end
    DESIGN.VOLS = num2str(tmp.dim(5));
    DESIGN.STANDARD = DESIGN.STRUCT;
    DESIGN.TOTAL_VOXELS = num2str(tmp.dim(2)*tmp.dim(3)*tmp.dim(4)*tmp.dim(5));
    DESIGN.FEAT_DIR = fullfile(session_dir,d{r},[func '.nii.gz']);
    DESIGN.CONDITION_1 = fullfile(session_dir,d{r},'condition_1.txt');
    DESIGN.CONDITION_2 = fullfile(session_dir,d{r},'condition_2.txt');
    DESIGN.CONDITION_3 = fullfile(session_dir,d{r},'condition_3.txt');
    DESIGN.CONDITION_4 = fullfile(session_dir,d{r},'condition_4.txt');
    DESIGN.CONDITION_5 = fullfile(session_dir,d{r},'condition_5.txt');
    DESIGN.CONDITION_6 = fullfile(session_dir,d{r},'condition_6.txt');
    DESIGN.CONDITION_7 = fullfile(session_dir,d{r},'task_trials.txt');
    fin = fopen(design_file,'rt');
    outname = [func '.feat_TTF.fsf'];
    fout = fopen(fullfile(session_dir,d{r},outname),'wt');
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
    commandc = ['feat ' fullfile(session_dir,d{r},outname)];
    disp(commandc);
    % Run feat
    if run_feat
        [~,~] = system([commandc ' &']);
    end
    % Also save this out in a file.
    fid = fopen(fullfile(session_dir, 'feat_TTF_scripts'), 'a');
    fprintf(fid, '%s\n',commandc);
    fclose(fid);
end
