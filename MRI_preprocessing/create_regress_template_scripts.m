function create_regress_template_scripts(session_dir,templateType,outDir,runs,func,saveDir,tcPart,leaveOut,V2V3)

% Creates shell scripts for running 'regress_template' on the UPenn
% cluster.
%
%   Usage:
%   create_regress_template_scripts(session_dir,templateType,outDir,runs,func,saveDir,tcPart,leaveOut)
%
%   Example:
%       session_dir = '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/';
%       templateType = 'coarse';
%       outDir = fullfile('/data/jet/abock/cluster_shell_scripts/fit_templates/ASB/10272014',templateType);
%       runs = '[2,4,6]'; % must be a string (!)
%       func = 's5.dbrf.tf';
%       tcPart = 'full';
%       leaveOut = '[1 3]';
%       V2V3 = '0'; % Don't use V2-V3 correlations
%       saveDir = fullfile(session_dir,'pRFs',templateType,func,'Movie');
%
%   Written by Andrew S Bock Sep 2015

%% Set defaults
hemis = {'lh' 'rh'};
cluster = '1'; % run on cluster
if ~exist(outDir,'dir')
    mkdir(outDir);
end
if ~exist('tcPart','var')
    tcPart = 'full'; % 'half' splits into 1st and 2nd halves
end
if ~exist('leaveOut','var')
    leaveOut = '0'; % only relevant when "tcPart = 'half'", dictates which run halves to leave out.
end
%% Create regress_template scripts
% psiParams = {...
%     -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ...
%     };
% FCxParams = {...
%     -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 ...
%     };
% FCyParams = {...
%     -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5 ...
%     };
for hh = 1:length(hemis);
    hemi = hemis{hh};
    switch templateType
        case 'coarse'
            if strcmp(hemi,'lh')
                ii = 4:11;% 0.0 - 0.7
                jj = 4:11; % -0.4 - 0.3
                kk = 3:10; % -0.6 - 0.1
                %                 ii = 6:12;%  0.2 <->  0.8
                %                 jj = 3:9; % -0.5 <->  0.1
                %                 kk = 3:9; % -0.6 <->  0.0
            else
                ii = 1:8; % -0.3 - 0.4
                jj = 3:10; % -0.5 - 0.2
                kk = 7:14; % -0.2 - 0.5
                %                 ii = 2:8; % -0.2 <->  0.4
                %                 jj = 2:8; % -0.6 <->  0.0
                %                 kk = 8:14;% -0.1 <->  0.5
            end
        case 'fine'
            ii = 1:5;
            jj = 1:5;
            kk = 1:5;
    end
    if strcmp(templateType,'pRF') || strcmp(templateType,'anat');
        job_name = [hemi '.' templateType '.regress'];
        matlab_string = ([...
            'regress_template(''' session_dir ''',''' saveDir ''',''' templateType ''',' ...
            runs ',''' hemi ''',''' func ''',''' templateType ''',' cluster ',''' tcPart ''',' leaveOut ',' V2V3 ');']);
        create_job_shell(outDir,job_name,matlab_string)
    else
        for i = ii
            for j = jj
                for k = kk
                    job_name = [hemi '.' num2str(i) '.' num2str(j) '.' num2str(k) '.regress'];
                    matlab_string = ([...
                        'regress_template(''' session_dir ''',''' saveDir ''',''' templateType ''',' ...
                        runs ',''' hemi ''',''' func ''',''' ...
                        [num2str(i) '.' num2str(j) '.' num2str(k)] ''',' ...
                        cluster ',''' tcPart ''',' leaveOut ',' V2V3 ');']);
                    create_job_shell(outDir,job_name,matlab_string)
                end
            end
        end
    end
end
%% Create submit script
for hh = 1:length(hemis);
    hemi = hemis{hh};
    submit_name = ['submit_' hemi '_regress'];
    job_string = listdir(fullfile(outDir,[hemi '*regress*.sh']),'files');
    cores = 1;
    mem = 5;
    create_submit_shell(outDir,submit_name,job_string,cores,mem);
end
system(['chmod +x ' fullfile(outDir,'*')]);