function create_regress_template_scripts(session_dir,templateType,outDir,runs,func,saveDir,tcPart,leaveOut,V2V3,logDir)

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
if ~exist('logDir','var')
    logDir = '/data/jet/abock/LOGS';
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
for hhh = 1:length(hemis);
    hemi = hemis{hhh};
    switch templateType
        case 'coarse'
            hh = 1:7;
            ii = 1:7;
            jj = 1:7;
        case 'fine'
            hh = 1:5;
            ii = 1:5;
            jj = 1:5;
        case 'coarseV2V3size'
            hh = 1:7;
            ii = 1:7;
            jj = 1:7;
            kk = 1:7;
            ll = 1:7;
    end
    if strcmp(templateType,'pRF') || strcmp(templateType,'anat');
        job_name = [hemi '.' templateType '.regress'];
        matlab_string = ([...
            'regress_template(''' session_dir ''',''' saveDir ''',''' templateType ''',' ...
            runs ',''' hemi ''',''' func ''',''' templateType ''',' cluster ',''' tcPart ''',' leaveOut ',' V2V3 ');']);
        create_job_shell(outDir,job_name,matlab_string);
    elseif strcmp(templateType,'coarseV2V3size')
        for h = hh
            for i = ii
                % conditional statement for the fact that large V2 AND V3 
                % templates could not be created (issue with model in Mathematica)
                if i < 6 || i == 6 && h < 7 
                    for j = jj
                        for k = kk
                            for l = ll
                                job_name = [hemi '.' num2str(h) '.' num2str(i) '.' ...
                                    num2str(j) '.' num2str(k) '.' num2str(l) '.regress'];
                                matlab_string = ([...
                                    'regress_template(''' session_dir ''',''' saveDir ''',''' templateType ''',' ...
                                    runs ',''' hemi ''',''' func ''',''' ...
                                    [num2str(h) '.' num2str(i) '.' num2str(j) '.' num2str(k) '.' num2str(l)] ''',' ...
                                    cluster ',''' tcPart ''',' leaveOut ',' V2V3 ');']);
                                create_job_shell(outDir,job_name,matlab_string);
                            end
                        end
                    end
                end
            end
        end
    else
        for h = hh
            for i = ii
                for j = jj
                    job_name = [hemi '.' num2str(h) '.' num2str(i) '.' ...
                        num2str(j) '.regress'];
                    matlab_string = ([...
                        'regress_template(''' session_dir ''',''' saveDir ''',''' templateType ''',' ...
                        runs ',''' hemi ''',''' func ''',''' ...
                        [num2str(h) '.' num2str(i) '.' num2str(j)] ''',' ...
                        cluster ',''' tcPart ''',' leaveOut ',' V2V3 ');']);
                    create_job_shell(outDir,job_name,matlab_string);
                end
            end
        end
    end
end
%% Create submit script
for hhh = 1:length(hemis);
    hemi = hemis{hhh};
    submit_name = ['submit_' hemi '_regress'];
    job_string = listdir(fullfile(outDir,[hemi '*regress*.sh']),'files');
    mem = 5;
    create_submit_shell(outDir,logDir,submit_name,job_string,mem);
end