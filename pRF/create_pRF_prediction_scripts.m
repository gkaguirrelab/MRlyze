function create_pRF_prediction_scripts(outDir,logDir,imFiles,paramFiles,outFiles,mem)

% Creates shell scripts for running 'createPredpRF'
%
%   Usage:
%   create_pRF_prediction_scripts(outDir,logDir,imFiles,paramFiles,outFiles,mem)
%
%   Inputs:
%   outDir          = output directory for shell scripts
%   logDir          = log directory for cluster output logs
%   imFiles         = cell matrix of pRF image files (full path)
%   paramFiles      = cell matrix of pRF param files (full path)
%   outFiles        = cell matrix of output .mat files (full path)
%   mem             = memory for cluster (default = 15)
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
    mem = 15;
end
%% Create run_pRF
for i = 1:length(imFiles)
    job_name = ['createPredpRF_' num2str(i)];
    matlab_string = (['createPredpRF(''' imFiles{i} ''',''' ...
        paramFiles{i} ''',''' outFiles{i} ''');']);
    create_job_shell(outDir,job_name,matlab_string)
end
%% Create submit script
submit_name = 'submit_createPredpRF';
job_string = listdir(fullfile(outDir,'*.sh'),'files');
create_submit_shell(outDir,logDir,submit_name,job_string,mem);