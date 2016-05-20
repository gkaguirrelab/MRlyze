function create_pRF_scripts(outDir,logDir,session_dir,subject_name,runNums,hemis,srcROIs,funcs)

% Creates shell scripts for running 'run_pRF'
%
%   Usage:
%   create_pRF_scripts(outDir,logDir,session_dir,subject_name,runNums,hemis,srcROIs,funcs)
%
%   Inputs:
%   outDir          = output directory for shell scripts
%   logDir          = log directory for cluster output logs
%   session_dir     = session directory
%   subject_name    = freesurfer subject name
%   runNums         = bold run numbers (double, not a string)
%   hemis           = hemispheres (e.g. hemis = {'lh' 'rh'};)
%   srcROIs         = source ROIs (e.g. srcROIs = {'cortex' 'volume'};)
%   funcs           = functional volume (e.g. funcs = 'wdrf.tf';)
%
%   Written by Andrew S Bock May 2016

%% Set defaults
if ~exist(outDir,'dir')
    mkdir(outDir);
end
if isempty(logDir)
    logDir = '/data/jet/abock/LOGS';
end
%% Create run_pRF
for rr = 1:length(runNums)
    runNum = runNums(rr);
    for hh = 1:length(hemis);
        hemi = hemis{hh};
        for ro = 1:length(srcROIs)
            srcROI = srcROIs{ro};
            for ff = 1:length(funcs)
                func = funcs{ff};
                job_name = [hemi '.' num2str(runNum) '.' srcROI '.' func '.run_pRF'];
                matlab_string = (['run_pRF(''' session_dir ''',''' subject_name ''',' ...
                    num2str(runNum) ',''' hemi ''',''' srcROI ''',''' func ''');']);
                create_job_shell(outDir,job_name,matlab_string)
            end
        end
    end
end
%% Create submit script
submit_name = ['submit_' subject_name '_run_pRF'];
job_string = listdir(fullfile(outDir,'*.sh'),'files');
mem = 40;
create_submit_shell(outDir,logDir,submit_name,job_string,mem);