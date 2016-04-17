function create_pRF_scripts(outDir,session_dir,subject_name,runNums,hemis,srcROIs)

% Creates shell scripts for running 'run_pRF' and 'average_pRF' on the
%   UPenn cluster.
%
%   Usage:
%   create_CF_scripts(outDir,session_dir,subject_name,runNums,templates,hemis,srcROIs,srcfunc,trgfunc,cond,V1only,cluster)
%
%   Written by Andrew S Bock Jan 2016


%%


run_pRF(session_dir,subject_name,runNum,hemi,srcROI)


average_pRF(session_dir,subject_name,runNums,srcROI);


%% Set defaults
if ~exist(outDir,'dir')
    mkdir(outDir);
end
%% Create run_pRF
for rr = 1:length(runNums)
    runNum = runNums(rr);
    for hh = 1:length(hemis);
        hemi = hemis{hh};
        job_name = [hemi '.' num2str(runNum) '.run_pRF'];
        matlab_string = (['run_pRF(''' session_dir ''',''' subject_name ''',' ...
            num2str(runNum) ',''' hemi ''',''' srcROI ''');']);
        create_job_shell(outDir,job_name,matlab_string)
        
    end
end
%% Create submit script (make_decimated_CF_predictions)
for tt = 1:length(templates)
    template = templates{tt};
    submit_name = ['submit_make_decimated_CF_predictions_' template];
    job_string = listdir(fullfile(outDir,['*' template '*make_decimated_CF_predictions.sh']),'files');
    cores = 1;
    mem = 20;
    create_submit_shell(outDir,submit_name,job_string,cores,mem);
end
system(['chmod +x ' fullfile(outDir,'/*')]);
%% Create find_best_CF
for rr = 1:length(runNums)
    runNum = runNums(rr);
    for hh = 1:length(hemis);
        hemi = hemis{hh};
        for ss = 1:length(srcROIs)
            srcROI = srcROIs{ss};
            for tt = 1:length(templates)
                template = templates{tt};
                job_name = [hemi '.' srcROI '.' template '.' num2str(runNum) '.find_best_CF'];
                matlab_string = (['find_best_CF(''' session_dir ''',' num2str(runNum) ...
                    ',''' hemi ''',''' srcROI ''',''' template ''',''' srcfunc ''',''' ...
                    trgfunc ''',''' cond ''',' num2str(V1only) ',' num2str(cluster) ');']);
                create_job_shell(outDir,job_name,matlab_string)
            end
        end
    end
end
%% Create submit script (find_best_CF)
for tt = 1:length(templates)
    template = templates{tt};
    submit_name = ['submit_find_best_CF_' template];
    job_string = listdir(fullfile(outDir,['*' template '*find_best_CF.sh']),'files');
    cores = 1;
    mem = 20;
    create_submit_shell(outDir,submit_name,job_string,cores,mem);
end
system(['chmod +x ' fullfile(outDir,'/*')]);