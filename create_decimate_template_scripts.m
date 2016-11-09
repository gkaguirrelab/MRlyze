function create_decimate_template_scripts(params)

% Creates shell scripts for running 'decimate_templates'
%
%   Usage:
%   create_decimate_template_scripts(params)
%
%   Inputs:
%   params.outDir           - output directory for shell scripts
%   params.logDir           - log directory for cluster output logs
%   params.subjectName      - freesurfer subject name
%   params.templateDir      - directory with the .mgz templates
%   params.mem              - memory for cluster (default = 10)
%
%   Written by Andrew S Bock Nov 2016

%% make the job script
job_name = [params.subjectName '.decimateTemplate'];
matlab_string = (['decimate_templates(''' params.subjectName ''',''' ...
    params.templateDir ''');']);
create_job_shell(params.outDir,job_name,matlab_string)
%% create the submit script
submit_name = ['submit.' params.subjectName '.decimateTemplate'];
job_string  = {fullfile(params.outDir,[job_name '.sh'])};
create_submit_shell(params.outDir,params.logDir,submit_name,job_string,params.mem);