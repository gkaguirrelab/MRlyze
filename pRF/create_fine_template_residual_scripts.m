function create_fine_template_residual_scripts(session_dir,outDir,runs,email)

% Creates shell scripts for running 'calc_fine_template_residual' on the UPenn
% cluster.
%
%   Usage:
%   create_fine_template_residual_scripts(session_dir,outDir,runs,email)
%
%   Example:
%   session_dir = '/data/jet/abock/data/Retinotopy/AEK/10012014/';
%   outDir = '/data/jet/abock/cluster_shell_scripts/calc_stimulus_corr/AEK/';
%   runs = '[2,4,6]'; % must be a string (!)
%
%   Written by Andrew S Bock Oct 2015

%% Set defaults
hemis = {'lh' 'rh'};
%% Create regress_template scripts
for hh = 1:length(hemis);
    hemi = hemis{hh};
    for i = 1:5
        for j = 1:5
            for k = 1:5
                job_name = [hemi '.' num2str(i) '.' num2str(j) '.' num2str(k) '.calc_fine_template_residual'];
                matlab_string = ([...
                    'calc_fine_template_residual(''' session_dir ''',' ...
                    runs ',''' hemi ''',''' [num2str(i) '.' num2str(j) '.' num2str(k)] ''');' ...
                    ]);
                create_job_shell(outDir,job_name,matlab_string)
            end
        end
    end
end
%% Create submit script
for hh = 1:length(hemis);
    hemi = hemis{hh};
    submit_name = ['submit_' hemi '_calc_fine_template_residual'];
    job_string = listdir(fullfile(outDir,[hemi '*calc_fine_template_residual*.sh']),'files');
    cores = 1;
    mem = 3;
    if ~exist('email','var');
        create_submit_shell(outDir,submit_name,job_string,cores,mem);
    else
        create_submit_shell(outDir,submit_name,job_string,cores,mem,email);
    end
end
cd(outDir);
system('chmod +x ./*');