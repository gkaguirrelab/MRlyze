function create_pRF_calc_scripts(outDir,logDir,outFiles,predFiles,inFiles,srcInds,mem)

% Creates shell scripts for running 'createPredpRF'
%
%   Usage:
%   create_pRF_calc_scripts(outDir,logDir,outFile,predFile,inFile,srcInds,mem)
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
%% Create run_pRF
for i = 1:length(outFiles)
    job_name = sprintf('calcpRF_%04d',i);
    matlab_string = (['calcpRF(''' outFiles{i} ''',''' ...
        predFiles{i} ''',''' inFiles{i} ''',[' (srcInds{i}) ']);']);
    create_job_shell(outDir,job_name,matlab_string)
end
%% Create submit script
submit_name = 'submit_calcpRF';
job_string = fullfile(outDir,listdir(fullfile(outDir,'*.sh'),'files'));
create_submit_shell(outDir,logDir,submit_name,job_string,mem);