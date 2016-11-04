function create_convert_Mathematica_template_scripts(outDir,logDir,session_dir,template_dir,mem)

% Creates shell scripts for running 'convert_Mathematica_template'
%
%   Usage:
%   create_pRF_scripts(outDir,logDir,session_dir,subject_name,runNums,hemis,srcROIs,funcs,mem)
%
%   Inputs:
%   outDir          = output directory for shell scripts
%   logDir          = log directory for cluster output logs
%   session_dir     = session directory
%   template_dir    = directory with the .mgz templates
%   mem             = memory for cluster (default = 10)
%
%   Written by Andrew S Bock May 2016

%% Set defaults
if ~exist(outDir,'dir')
    mkdir(outDir);
end
if isempty(logDir)
    logDir = '/data/jet/abock/LOGS';
end
if ~exist('mem','var')
    mem = 10;
end
%% Find template files
t = listdir(fullfile(template_dir,'*.mgz'),'files');
%% Create run_pRF
for tt = 1:length(t)
    job_name = sprintf('convertTemplate.%05d',tt);
    matlab_string = (['convert_Mathematica_template(''' session_dir ''',''' ...
        template_dir ''',' num2str(tt) ');']);
    create_job_shell(outDir,job_name,matlab_string)
end
%% Create submit script
submit_name = 'submit_convertTemplate';
job_string = listdir(fullfile(outDir,'*.sh'),'files');
create_submit_shell(outDir,logDir,submit_name,job_string,mem);