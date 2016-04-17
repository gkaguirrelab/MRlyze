function feat_higher_level_4feat_5copes(session_dir,subject_name,runNums,func,run_feat,design_file)

% Creates a 'feat_higher_level_4feat_5copes.fsf' file in each bold directory for a given
% session directory.
%
%   Usage:
%   feat_higher_level_4feat_5copes(session_dir,subject_name,runNums,func,run_feat,design_file)
%
%   Defaults:
%     func = 's2.drf.tf'; % default functional input
%     run_feat = 0; % default will NOT run feat
%     design_file = 'feat_higher_level_4feat_5copes_template.fsf';
%
%   Written by Andrew S Bock Oct 2015

%% Set default parameters
if ~exist('func','var')
    func = 's2.drf.tf';
end
if ~exist('run_feat','var')
    run_feat = 0; % run FEAT from matlab.  Alternative is to run from the terminal, e.g. if on a cluster
end
if ~exist('design_file','var')
    design_dir = which('feat_higher_level_4feat_5copes');
    design_dir = fileparts(design_dir);
    design_file = fullfile(design_dir,'feat_higher_level_4feat_5copes_template.fsf');
end
%% Add to log
SaveLogInfo(session_dir,mfilename,session_dir,runNums,func,run_feat,design_file);
%% Find bold run directories
d = find_bold(session_dir);
%% Overwrite feat reg
overwrite_feat_reg(session_dir,subject_name,runNums,func);
%% Run FEAT
DESIGN.FEAT1 = fullfile(session_dir,d{runNums(1)},[func '.feat']);
DESIGN.FEAT2 = fullfile(session_dir,d{runNums(2)},[func '.feat']);
DESIGN.FEAT3 = fullfile(session_dir,d{runNums(3)},[func '.feat']);
DESIGN.FEAT4 = fullfile(session_dir,d{runNums(4)},[func '.feat']);
fin = fopen(design_file,'rt');
outname = [func '.feat_higher_level_4feat_5copes.fsf'];
fout = fopen(fullfile(session_dir,d{runNums(1)},outname),'wt');
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
disp('Feat commmand:')
commandc = ['feat ' fullfile(session_dir,d{runNums(1)},outname)];
disp(commandc);
% Run feat
if run_feat
    disp('Running above feat command...');
    [~,~] = system([commandc ' &']);
end
% Also save this out in a file.
fid = fopen(fullfile(session_dir, 'feat_higher_level_scripts'), 'a');
fprintf(fid, '%s\n',commandc);
fclose(fid);