%%
Control = {'A111907G';'D030208S';'L030208D';'R030308W';'S102907D';
    'W021808H';'M042507D'; 'R042507M';'S042507C';'S042507H'};
Blind = {'C111507D';'C111907L';'D010908G';'E011108K';'E122007P';'M012108K';
    'M032408K';'M110707N';'V020808H';'V061908W';'V020408W'};
subjects = [Control; Blind];
datadir = '/jet/aguirre/abock/Semantic_Decoding';
savedir = '~/data/Semantic_Decoding'; % can't save to /jet/aguirre
% 'searchlight_results_2_conditions_aud_tac'
% 'searchlight_results_2_conditions_aud_vis'
% 'searchlight_results_2_conditions_tac_vis'
% 'searchlight_results_3_conditions'
% 'searchlight_results_4_conditions'
comp_file = 'searchlight_results_2_conditions_tac_vis';
hemi = {'lh' 'rh'};
smooth_size = 10;
%% Confirm data exists for every subject
ct = 0;
for s = 1:length(subjects)
    cd(fullfile(datadir,subjects{s}));
    if exist(fullfile('./',[comp_file '.nii.gz']))
        disp(cd)
        ct = ct + 1;
    end
end
disp(num2str(ct));
%% Project to surface
progBar = ProgressBar(length(subjects),'Projecting to surface');
for s = 1:length(subjects)
    for h = 1:length(hemi)
        cd(fullfile(datadir,subjects{s}));
        [~,~] = system(['mri_vol2surf --mov ./' comp_file '.nii.gz --reg ./bbreg.dat' ...
            ' --hemi ' hemi{h} ' --projfrac 0.5 --surf-fwhm ' num2str(smooth_size) ...
            ' --o ./' hemi{h} '_' comp_file '_surf.nii.gz']);
        [~,~] = system(['mri_surf2surf --srcsubject ' subjects{s} ' --sval ./' ...
            hemi{h} '_' comp_file '_surf.nii.gz --trgsubject fsaverage_sym --hemi ' hemi{h} ...
            ' --tval ./' hemi{h} '_'  comp_file '_fssym_surf.nii.gz']);
    end
    progBar(s);
end
%% Controls
for h = 1:length(hemi)
    avgnii = [];
    for s = 1:length(Control)
        cd(fullfile(datadir,Control{s}));
        nii = load_nifti(['./' hemi{h} '_' comp_file '_fssym_surf.nii.gz']);
        avgnii = [avgnii nii.vol];
    end
    nii.vol = mean(avgnii,2);
    % Save Average Matrix
    save(fullfile(savedir,[hemi{h} '_avg_Control_' comp_file]),'avgnii')
    save_nifti(nii,fullfile(savedir,[hemi{h} '_avg_Control_' comp_file '.nii.gz']));
end
%% Blind
for h = 1:length(hemi)
    avgnii = [];
    for s = 1:length(Blind)
        cd(fullfile(datadir,Blind{s}));
        nii = load_nifti(['./' hemi{h} '_' comp_file '_fssym_surf.nii.gz']);
        avgnii = [avgnii nii.vol];
    end
    nii.vol = mean(avgnii,2);
    % Save Average Matrix
    save(fullfile(savedir,[hemi{h} '_avg_Blind_' comp_file]),'avgnii')
    save_nifti(nii,fullfile(savedir,[hemi{h} '_avg_Blind_' comp_file '.nii.gz']));
end
%% T-test
%% Controls
for h = 1:length(hemi)
    tnii = [];
    for s = 1:length(Control)
        cd(fullfile(datadir,Control{s}));
        nii = load_nifti(['./' hemi{h} '_' comp_file '_fssym_surf.nii.gz']);
        tnii = [tnii nii.vol];
    end
    tstat = nan(length(tnii),1);
    H = tstat;
    progBar = ProgressBar(length(tnii),'Calculating t-stats');
    for t = 1:length(tnii)
        clear stats
        if strcmp(comp_file,'searchlight_results_3_conditions');
            [H(t),~,~,stats] = ttest(tnii(t,:) - 1/3);
        elseif strcmp(comp_file,'searchlight_results_4_conditions')
            [H(t),~,~,stats] = ttest(tnii(t,:) - 1/4);
        end
        tstat(t) = stats.tstat;
        if ~mod(t,100);progBar(t);end
    end
    nii.vol = tstat;
    save_nifti(nii,fullfile(savedir,[hemi{h} '_tstat_Control_' comp_file '.nii.gz']));
    nii.vol = H;
    save_nifti(nii,fullfile(savedir,[hemi{h} '_H_Control_' comp_file '.nii.gz']));
end
%% Blind
for h = 1:length(hemi)
    tnii = [];
    for s = 1:length(Blind)
        cd(fullfile(datadir,Blind{s}));
        nii = load_nifti(['./' hemi{h} '_' comp_file '_fssym_surf.nii.gz']);
        tnii = [tnii nii.vol];
    end
    tstat = nan(length(tnii),1);
    H = tstat;
    progBar = ProgressBar(length(tnii),'Calculating t-stats');
    for t = 1:length(tnii)
        clear stats
        if strcmp(comp_file,'searchlight_results_3_conditions');
            [H(t),~,~,stats] = ttest(tnii(t,:) - 1/3);
        elseif strcmp(comp_file,'searchlight_results_4_conditions')
            [H(t),~,~,stats] = ttest(tnii(t,:) - 1/4);
        end
        tstat(t) = stats.tstat;
        if ~mod(t,100);progBar(t);end
    end
    nii.vol = tstat;
    save_nifti(nii,fullfile(savedir,[hemi{h} '_tstat_Blind_' comp_file '.nii.gz']));
    nii.vol = H;
    save_nifti(nii,fullfile(savedir,[hemi{h} '_H_Blind_' comp_file '.nii.gz']));
end
%% T-test LO
h = 1; % left hemi only
avgControl = load(fullfile(savedir,[hemi{h} '_avg_Control_' comp_file]),'avgnii');
avgBlind = load(fullfile(savedir,[hemi{h} '_avg_Blind_' comp_file]),'avgnii');
LO = load_nifti('~/data/LH-Visuotopic-LO.nii.gz');
LOind = find(LO.vol);
avgControl.LO = mean(avgControl.avgnii(LOind,:));
avgBlind.LO = mean(avgBlind.avgnii(LOind,:));
x = avgControl.LO';
y = avgBlind.LO';
mean(x)
mean(y)
[H,P,CI,stats] = ttest2(avgControl.LO',avgBlind.LO')




