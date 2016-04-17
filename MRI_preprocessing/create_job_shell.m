function create_job_shell(outDir,job_name,matlab_string)

% Writes out a shell script to run matlab code on the cluster.
%
%   Usage:
%   create_job_shell(outDir,job_name,matlab_string)
%
%   example:
%   outDir = '~/cluster_scripts';
%   job_name = 'job_example';
%   matlab_string = 'sort_nifti(''~/data/session_dir'');'
%   create_job_shell(outDir,job_name,matlab_string);
%
%   The job_name file ('e.g. job_example.sh') contains the matlab script
%   specified by 'matlab_script'.  In our example, this file looks like:
%
%   #!/bin/bash
%   matlab -nodisplay -nosplash -r "sort_nifti('~/data/session_dir');"
%
%   Written by Andrew S Bock Aug 2015

%% create shell script to run matlab code
fid = fopen(fullfile(outDir,[job_name '.sh']),'w');
%fprintf(fid,'#!/bin/bash\n');
full_matlab_string = (['matlab -nodisplay -r "' matlab_string '"']);
fprintf(fid,full_matlab_string);
fclose(fid);