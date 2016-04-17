function feat_higher_level(session_dirs,subject_name,dirs,func,design_file,SUBJECTS_DIR)

%   Following feat_stats, this function will run higher level stats across
%   the various runs/sessions, using the resulting feat directories from
%   feat_stats.
%
%   Usage:
%   feat_higher_level(session_dirs,subject_name,dirs,func,design_file,SUBJECTS_DIR)
%
%   inputs:
%   session_dirs - cell containing strings of the paths to the session
%   directories of interest
%   subject_name - freesurfer subject name
%   dirs - vector of feat directories (e.g. [1 4 7 10 13 16])
%   directories of interest
%   func - functional data file (default - 'dbrf.tf')
%
%   Written by Andrew S Bock Dec 2014
%
%
%   03/09/2015  ASB  created 'overwrite_feat_reg'

%% Set default parameters
if ~exist('session_dirs','var')
    error('"session_dirs" not defined')
end
if ~exist('subject_name','var')
    error('"subject_name" not defined')
end
if ~exist('dirs','var')
    error('"dirs" not defined') % feat directories
end
if ~exist('func','var')
    func = 'sdbrf.tf'; % functional data file
end
design_dir = which('feat_higher_level');
design_dir = fileparts(design_dir);
if ~exist('design_file','var')
    design_file = fullfile(design_dir,'feat_higher_level_template.fsf');
end
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
feat_dir = [func '.feat'];
ct = 0;
%% Add to log
SaveLogInfo(session_dirs{1},mfilename,cell2str(session_dirs),subject_name,dirs,func,design_file,SUBJECTS_DIR)

%% Overwrite feat registrations with bbregister registration
outdir = overwrite_feat_reg(session_dirs,subject_name,dirs,func);

%% Run FEAT
% Create 'dirs' variable, which is a list of all bold directories in the
% session_dirs variable
for s = 1:length(session_dirs)
    session_dir = session_dirs{s};
    % Find bold run directories
    d = listdir(fullfile(session_dir,'*bold_*'),'dirs');
    if isempty(d)
        d = listdir(fullfile(session_dir,'BOLD_*'),'dirs');
    end
    nruns = length(d);
    for r = 1:nruns
        ct = ct+1;
        feat_dirs{ct} = fullfile(session_dir,d{r},feat_dir);
    end
end
% Specify the specific directories in the 'feat_dirs' variable
for f=1:length(dirs)
    eval(['DESIGN.FEAT_' num2str(f) ' = feat_dirs{dirs(' num2str(f) ')};']);
end
% Fill out 'feat_higher_level.fsf'
fin = fopen(design_file,'rt');
fout = fopen(fullfile(outdir,'feat_higher_level.fsf'),'wt');
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
commandc = ['feat ' fullfile(outdir,'feat_higher_level.fsf')];
disp(commandc);
% Run feat
[~,~] = system([commandc ' &']);
% Also save this out in a file.
fid = fopen(fullfile(outdir, 'feat_higher_level_scripts'), 'a');
fprintf(fid, '%s\n',commandc);
fclose(fid);