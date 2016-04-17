function create_submit_shell(outDir,submit_name,job_string,cores,mem,email)

% Writes out a shell script to run matlab code on the cluster.
%
%   Usage:
%   create_submit_shell(outDir,submit_name,job_string,cores,mem,email)
%
%   example:
%   outDir = '~/cluster_scripts';
%   submit_name = 'submit_example'
%   job_string = 'job_example_1.sh'; % this is a cell, and can contain many
%       job strings
%   cores = 4;
%   mem = 16;
%   create_submit_shell(outDir,submit_name,job_string,cores,mem)
%
%   This creates a file 'submit_example.sh' in outDir, with one line:
%       $SGE_ROOT/bin/linux-x64/qsub -binding linear:4 -pe unihost 4 -l
%       h_vmem=16.2G,s_vmem=16G job_example_1.sh
%
%   Written by Andrew S Bock Aug 2015

%% Set defaults
if ~exist('cores','var');
    cores = num2str(1);
end
if ~exist('mem','var');
    mem = 10;
end
h_vmem = ([num2str(mem) '.2G']);
s_vmem = ([num2str(mem) 'G']);

%% create shell string to submit jobs
fid = fopen(fullfile(outDir,[submit_name '.sh']),'w');
for L = 1:length(job_string);
    if ~exist('email','var');
        shell_string = (['$SGE_ROOT/bin/lx-amd64/qsub -e /data/jet/abock/LOGS ' ...
            '-o /data/jet/abock/LOGS -binding linear:' num2str(cores) ...
            ' -pe unihost ' num2str(cores) ' -l h_vmem=' h_vmem ',s_vmem=' s_vmem ...
            ' ' job_string{L}]);
    else
        shell_string = (['$SGE_ROOT/bin/lx-amd64/qsub -e /data/jet/abock/LOGS ' ...
            '-o /data/jet/abock/LOGS -binding linear:' num2str(cores) ...
            ' -pe unihost ' num2str(cores) ' -l h_vmem=' h_vmem ',s_vmem=' s_vmem ...
            ' -m e -M ' email ' ' job_string{L}]);
    end
    %     if ~exist('email','var');
    %     shell_string = (['qsub -l h_vmem=' h_vmem ',s_vmem=' s_vmem ...
    %         ' ' job_string{L}]);
    %     else
    %         shell_string = (['qsub -l h_vmem=' h_vmem ',s_vmem=' s_vmem ...
    %         ' -m e -M ' email ' ' job_string{L}]);
    %     end
    fprintf(fid,[shell_string '\n']);
end
fclose(fid);