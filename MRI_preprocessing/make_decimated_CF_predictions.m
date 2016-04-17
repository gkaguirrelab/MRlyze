function [shiftPredResp,tc_shift,xList,sigList,V1trgdists,V2trgdists,V3trgdists] = ...
    make_decimated_CF_predictions(session_dir,subject_name,runNum,hemi,template,trgfunc,cond,DoG)

%% Creates connective field (CF) predictions
%   Similar to the pRF method, we create a set of CF predicted timecourses,
%   and will later search through this space (using find_best_CF.m)
%
%   Usage:
%   [shiftPredResp,tc_shift,xList,sigList,V1trgdists,V2trgdists,V3trgdists] = ...
%    make_decimated_CF_predictions(session_dir,subject_name,runNum,hemi,template,trgfunc,cond,DoG)
%
%   Written by Andrew S Bock Jan 2016
%% Set up defaults
if ~exist('hemi','var')
    hemi = 'lh';
end
if ~exist('template','var')
    template='fine';
end
if ~exist('trgfunc','var')
    trgfunc = 's5.dbrf.tf';
end
if ~exist('cond','var')
    cond = 'Movie';
end
if ~exist('DoG','var');
    DoG = 1;
end
decimation_level='0.1';
src_surf = 'inflated';
trg_surf = [decimation_level '.' src_surf];
%% Setup sigmas
% min/lin_step/max sigma in degree visual angle
ctrSigbnd = [1/4 1/4 20];
% Linear scale for visual field less than 5 degrees
Sigborder = 5;
% number of log sigmas
nLogSigs = 6;
% scale of surround sigma
srdSigma = 0; % (1.5:.75:3)';
% center and surround amplitudes
ctrBetas = 1; % (0:.25:1)';
srdBetas = 0; % (0:.25:1)';
%% Make sigma list
tmps = ctrSigbnd(1):ctrSigbnd(2):Sigborder;
linsigs = tmps(1:end-1);
logsigs = logspace(log10(Sigborder),log10(ctrSigbnd(3)),nLogSigs);
ctrSigma = [linsigs logsigs];
%% Find bold run directories
d = find_bold(session_dir);
nruns = length(d);
disp(['Session_dir = ' session_dir]);
disp(['Found ' num2str(nruns) ' runs']);
disp(['Working on ' d{runNum}]);
%% Decimate pRF templates
switch template
    case 'V1'
        % need to add this, currently don't have a decimated directory
    case 'pRF'
        areas = load_nifti(fullfile(session_dir,'pRFs','pRF_templates',...
            'decimated_templates',[hemi '.areas.pRF.nii.gz']));
        ecc = load_nifti(fullfile(session_dir,'pRFs','pRF_templates',...
            'decimated_templates',[hemi '.ecc.pRF.nii.gz']));
        pol = load_nifti(fullfile(session_dir,'pRFs','pRF_templates',...
            'decimated_templates',[hemi '.pol.pRF.nii.gz']));
    case 'fine'
        tdir = fullfile(session_dir,'pRFs','fine',trgfunc,cond);
        [~,~,sorted_templates] = find_best_template(template,tdir,hemi);
        bestTemplate = sorted_templates{1}(4:(strfind(sorted_templates{1},'varexp')-2));
        areas = load_nifti(fullfile(session_dir,'pRFs','fine_model_templates',...
            'decimated_templates',[hemi '.areas.' bestTemplate '.nii.gz']));
        ecc = load_nifti(fullfile(session_dir,'pRFs','fine_model_templates',...
            'decimated_templates',[hemi '.ecc.' bestTemplate '.nii.gz']));
        pol = load_nifti(fullfile(session_dir,'pRFs','fine_model_templates',...
            'decimated_templates',[hemi '.pol.' bestTemplate '.nii.gz']));
end
% Get target indices, ecc, pol for V1
V1trgind = find(abs(areas.vol)==1);
V1trgecc = ecc.vol(V1trgind);
V1trgpol = pol.vol(V1trgind);
% Get target indices, ecc, pol for V2
V2trgind = find(abs(areas.vol)==2);
V2trgecc = ecc.vol(V2trgind);
V2trgpol = pol.vol(V2trgind);
% Get target indices, ecc, pol for V3
V3trgind = find(abs(areas.vol)==3);
V3trgecc = ecc.vol(V3trgind);
V3trgpol = pol.vol(V3trgind);
%% Project sagittal sinus mask
% invol = fullfile(session_dir,'sagittal_sinus_mask.nii.gz');
% targvol = fullfile(session_dir,d{runNum},'single_TR.nii.gz');
% outvol = fullfile(session_dir,d{runNum},'sagittal_sinus_mask.nii.gz');
% reg = fullfile(session_dir,d{runNum},'brf_bbreg.dat');
% interp = 'nearest';
% inv = 1;
% mri_vol2vol(targvol,invol,outvol,reg,interp,inv);
%% Decimate target file
in_vol = fullfile(session_dir,d{runNum},[trgfunc '.surf.' hemi '.nii.gz']);
out_vol = fullfile(session_dir,d{runNum},[trgfunc '.surf.' trg_surf '.' hemi '.nii.gz']);
decimate_surf(subject_name,hemi,src_surf,trg_surf,in_vol,out_vol);
% Load target surface file
trgfile = out_vol;
%% Load timecourses
disp('Loading timecourses...');
% load files
trg = load_nifti(trgfile);
trgdims = size(trg.vol);
trgtc = reshape(trg.vol,trgdims(1)*trgdims(2)*trgdims(3),trgdims(4))';
% Pull out relevant timecourses (V1)
V1trgtc = trgtc(:,V1trgind);
V1trgtc = V1trgtc - repmat(mean(V1trgtc),size(V1trgtc,1),1);
% Pull out relevant timecourses (V2)
V2trgtc = trgtc(:,V2trgind);
V2trgtc = V2trgtc - repmat(mean(V2trgtc),size(V2trgtc,1),1);
% Pull out relevant timecourses (V3)
V3trgtc = trgtc(:,V3trgind);
V3trgtc = V3trgtc - repmat(mean(V3trgtc),size(V3trgtc,1),1);
% Set timecourses with very little variation (var<.1) to flat
V1trgtc = set_to_flat(V1trgtc);
V2trgtc = set_to_flat(V2trgtc);
V3trgtc = set_to_flat(V3trgtc);
disp('done.');
%% Calculate the average V1-V3 signal (include time shift)
% V1_V3tmpPredResp = mean([V1trgtc,V2trgtc,V3trgtc],2);
% V1_V3tmpPredResp = V1_V3tmpPredResp - mean(V1_V3tmpPredResp);
% %% Get indices and timecourse for the sagittal sinus
% sag_mask = load_nifti(fullfile(session_dir,d{runNum},'sagittal_sinus_mask.nii.gz'));
% sag_vol = load_nifti(fullfile(session_dir,d{runNum},[func '.nii.gz']));
% sagdims = size(sag_vol.vol);
% sag_tmptc = reshape(sag_vol.vol,sagdims(1)*sagdims(2)*sagdims(3),sagdims(4))';
% sagtc = sag_tmptc(:,find(sag_mask.vol==1));
% sagtc = sagtc - repmat(mean(sagtc),size(sagtc,1),1);
% sagSinustmpPredResp = mean(sagtc,2);
% sagSinustmpPredResp = sagSinustmpPredResp - mean(sagSinustmpPredResp);
% %% Regress out average V1-V3 and sagittal sinus from
% [regmat] = [ones(size(V1trgtc,1),1),V1_V3tmpPredResp,sagSinustmpPredResp];
% % V1
% for i = 1:size(V1trgtc,2);
%     [rbeta] = regress(V1trgtc(:,i),regmat);
%     V1trgtc(:,i) = V1trgtc(:,i)-regmat*(rbeta);
%     V1trgtc(:,i) = V1trgtc(:,i) - mean(V1trgtc(:,i));
% end
% % V2
% for i = 1:size(V2trgtc,2);
%     [rbeta] = regress(V2trgtc(:,i),regmat);
%     V2trgtc(:,i) = V2trgtc(:,i)-regmat*(rbeta);
%     V2trgtc(:,i) = V2trgtc(:,i) - mean(V2trgtc(:,i));
% end
% % V3
% for i = 1:size(V3trgtc,2);
%     [rbeta] = regress(V3trgtc(:,i),regmat);
%     V3trgtc(:,i) = V3trgtc(:,i)-regmat*(rbeta);
%     V3trgtc(:,i) = V3trgtc(:,i) - mean(V3trgtc(:,i));
% end
%% Set the tc_shift
% postive values move predicted responses earlier in time, negative values
% move later in time
TR = trg.pixdim(5)/1000;
disp(['TR = ' num2str(TR) ' seconds']);
tc_shift = (-1:0.5:5)/TR;
%% Get distances in target ROI in visual space (from V1 centers)
% V1
[V1trgX,V1trgY] = pol2cart(V1trgpol,V1trgecc);
for i = 1:length(V1trgind)
    V1trgdists(i,:) = sqrt( (V1trgX(i) - V1trgX).^2 + (V1trgY(i) - V1trgY).^2);
end
% V2
[V2trgX,V2trgY] = pol2cart(V2trgpol,V2trgecc);
for i = 1:length(V2trgind)
    V2trgdists(i,:) = sqrt( (V2trgX(i) - V1trgX).^2 + (V2trgY(i) - V1trgY).^2);
end
% V3
[V3trgX,V3trgY] = pol2cart(V3trgpol,V3trgecc);
for i = 1:length(V3trgind)
    V3trgdists(i,:) = sqrt( (V3trgX(i) - V1trgX).^2 + (V3trgY(i) - V1trgY).^2);
end
%% Create full sigma list
% Create sigma list for each center location
clear sigList
seedX = 1:length(V1trgdists);
[xList] = ndgrid(seedX,ctrSigma,srdSigma,ctrBetas,srdBetas);
xList = xList(:);
ct = 0;
for s1 = 1:length(ctrSigma);
    for s2 = 1:length(srdSigma);
        for s3 = 1:length(ctrBetas);
            for s4 = 1:length(srdBetas);
                ct = ct+1;
                singleSigList(ct,:) = [ctrSigma(s1),ctrSigma(s1)*srdSigma(s2),ctrBetas(s3),srdBetas(s4)];
            end
        end
    end
end
% Replicate this sigma list for all center locations
sigList = repmat(singleSigList,length(xList)/size(singleSigList,1),1);
% Make single (to save memory)
sigList = single(sigList);
xList = single(xList);
% Remove flat sigmas
badind = sigList(:,3) == 0 & sigList(:,4) == 0;
sigList(badind,:) = [];
xList(badind) = [];
%% Break up sigma list matrix
[sigPerTask,sigTasks] = calc_tasks(length(sigList),ceil(length(sigList)/1000));
sigidx = [];
for i = 1:sigTasks
    if isempty(sigidx);
        sigidx = [1,sigPerTask(i)];
    else
        sigidx = [sigidx;[sigidx(end,2)+1,sigidx(end,2)+sigPerTask(i)]];
    end
    sigvals{i} = sigidx(i,1):sigidx(i,2);
end
%% Make predicted timecourses
shiftPredResp = single(nan(3,length(tc_shift),size(V1trgtc,1),length(sigList)));
avgV1_V3tc = single(nan(length(tc_shift),size(V1trgtc,1)));
sagSinustc = avgV1_V3tc;
progBar = ProgressBar(length(tc_shift),'Making predictions...');
for tt = 1:length(tc_shift)
    for ss = 1:sigTasks;
        % Get centers
        V1x = V1trgdists(:,xList(sigvals{ss}));
        V2x = V2trgdists(:,xList(sigvals{ss}));
        V3x = V3trgdists(:,xList(sigvals{ss}));
        % Get sigmas (visual angle)
        V1sig = sigList(sigvals{ss},:);
        V2sig = sigList(sigvals{ss},:);
        V3sig = sigList(sigvals{ss},:);
        % Make predictions for V1
        [~,V1tmpPredResp] = calc_CF_pred(V1x,V1sig,V1trgtc,[],DoG);
        s = 1:size(V1tmpPredResp,1); % get the time
        u = s + tc_shift(tt); % shift the time (positive = shift left)
        tmp1 = sinc_interp(V1tmpPredResp',s,u); % sinc interpolate to shifted time
        tmp1 = tmp1'; % fix transpose from 'sinc_interp';
        shiftPredResp(1,tt,:,sigvals{ss}) = tmp1 - repmat(mean(tmp1),size(tmp1,1),1);
        % Make predictions for V2 (with closest V1 visual field center)
        [~,V2tmpPredResp] = calc_CF_pred(V2x,V2sig,V2trgtc,[],DoG);
        s = 1:size(V2tmpPredResp,1); % get the time
        u = s + tc_shift(tt); % shift the time (positive = shift left)
        tmp2 = sinc_interp(V2tmpPredResp',s,u); % sinc interpolate to shifted time
        tmp2 = tmp2'; % fix transpose from 'sinc_interp';
        shiftPredResp(2,tt,:,sigvals{ss}) = tmp2 - repmat(mean(tmp2),size(tmp2,1),1);
        % Make predictions for V3 (with closest V1 visual field center)
        [~,V3tmpPredResp] = calc_CF_pred(V3x,V3sig,V3trgtc,[],DoG);
        s = 1:size(V3tmpPredResp,1); % get the time
        u = s + tc_shift(tt); % shift the time (positive = shift left)
        tmp3 = sinc_interp(V3tmpPredResp',s,u); % sinc interpolate to shifted time
        tmp3 = tmp3'; % fix transpose from 'sinc_interp';
        shiftPredResp(3,tt,:,sigvals{ss}) = tmp3 - repmat(mean(tmp3),size(tmp3,1),1);
    end
    progBar(tt);
end
%% Save file
disp(['Saving file - ' fullfile(session_dir,d{runNum},...
    ['CF_predictions.' hemi '.' template '.' trgfunc '.' cond '.mat'])]);
save(fullfile(session_dir,d{runNum},...
    ['CF_predictions.' hemi '.' template '.' trgfunc '.' cond '.mat']),...
    'shiftPredResp','tc_shift','xList','sigList','V1trgdists','V2trgdists','V3trgdists','-v7.3');
disp('done.');