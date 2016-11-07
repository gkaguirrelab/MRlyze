function create_submit_shell(outDir,logDir,submit_name,job_string,mem)

% Writes out a shell script to run matlab code on the cluster.
%
%   Usage:
%   create_submit_shell(outDir,logDir,submit_name,job_string,mem)
%
%   example:
%   outDir = '/some/directory/for/scripts/'
%   logDir = '/data/jet/abock/LOGS';
%   submit_name = 'submit_example'
%   job_string = 'job_example_1.sh'; % this is a cell, and can contain many
%       job strings
%   mem = 5;
%   create_submit_shell(outDir,logDir,submit_name,job_string,mem)
%
%   Written by Andrew S Bock Aug 2015

%% Set defaults
if ~exist('mem','var');
    mem = 5;
end
h_vmem = ([num2str(mem) '.2G']);
s_vmem = ([num2str(mem) 'G']);

%% create shell string to submit jobs
fid = fopen(fullfile(outDir,[submit_name '.sh']),'w');
for L = 1:length(job_string);
    shell_string = (['qsub -e ' logDir ...
        ' -o ' logDir ' -l h_vmem=' h_vmem ',s_vmem=' s_vmem ...
        ' ' job_string{L}]);
    fprintf(fid,[shell_string '\n']);
end
fclose(fid);