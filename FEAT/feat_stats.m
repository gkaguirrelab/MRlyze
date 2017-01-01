function feat_stats(session_dir,rr,func,standard,design_dir,design_file,run_feat)

%   Following denoise, this function will run stats and poststats in
%   feat using the nuisance regressors from motion and physio, as well as
%   GLM using custom text files for each condition.
%
%   Written by Andrew S Bock Oct 2014

%% Set default parameters
if ~exist('session_dir','var')
    error('"session_dir" not defined')
end
if ~exist('rr','var')
    error('"run" not defined')
end
if ~exist('func','var')
    func = 'sdbrf.tf'; % functional data file
end
if ~exist('standard','var')
    standard = fullfile(session_dir,'MPRAGE','001','MPRAGE_brain.nii.gz'); % functional data file
end
if ~exist('design_dir','var')
    design_dir = which('feat_stats');
    design_dir = fileparts(design_dir);
end
if ~exist('design_file','var')
    design_file = fullfile(design_dir,'feat_stats_template.fsf');
end
if ~exist('run_feat','var')
    run_feat = 1; % default to also run feat
end
%% Add to log
SaveLogInfo(session_dir,mfilename,session_dir,rr,func,standard,design_dir,design_file)

%% Find bold run directories
d = listdir(fullfile(session_dir,'*bold_*'),'dirs');
if isempty(d)
    d = listdir(fullfile(session_dir,'BOLD_*'),'dirs');
end
%% Run FEAT
disp('Loading functional file...');
tmp = load_nifti(fullfile(session_dir,d{rr},[func '.nii.gz']));
disp('done.');
DESIGN.TR = num2str(tmp.pixdim(5)/1000); % TR is in msec, convert to sec
if tmp.pixdim(5) < 100 % use 100, in case very short TR is used (i.e. multi-band)
    error('TR is not in msec');
end
DESIGN.VOLS = num2str(tmp.dim(5));
DESIGN.STANDARD = standard;
DESIGN.TOTAL_VOXELS = num2str(tmp.dim(2)*tmp.dim(3)*tmp.dim(4)*tmp.dim(5));
DESIGN.FEAT_DIR = fullfile(session_dir,d{rr},func);
% Loop through condition text files
c = listdir(fullfile(session_dir,d{rr},'condition_*.txt'),'files');
for cc = 1:length(c)
    eval(['DESIGN.CONDITION_' num2str(cc) ' = fullfile(session_dir,' ...
        'd{rr},''condition_' num2str(cc) '.txt'');']);
end
% Write out task trials
eval(['DESIGN.CONDITION_' num2str((length(c)+1)) ' = fullfile(session_dir,' ...
    'd{rr},''task_trials.txt'');']);
% Find if there were any bad trials
badTrials = fullfile(session_dir,d{rr},'bad_trials.txt');
if exist(badTrials,'file');
    eval(['DESIGN.CONDITION_' num2str((length(c)+2)) ' = fullfile(session_dir,' ...
        'd{rr},''bad_trials.txt'');']);
end
% Open and write to design file
fin = fopen(design_file,'rt');
fout = fopen(fullfile(session_dir,d{rr},'feat_stats.fsf'),'wt');
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
commandc = ['feat ' fullfile(session_dir,d{rr},'feat_stats.fsf')];
disp(commandc);
% Run feat
if run_feat
    [~,~] = system([commandc ' &']);
end
% Also save this out in a file.
fid = fopen(fullfile(session_dir, 'feat_stats_scripts'), 'a');
fprintf(fid, '%s\n',commandc);
fclose(fid);