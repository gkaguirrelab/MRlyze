function [shiftPredResp,tc_shift,xList,sigList,V1trgdists,V2trgdists,V3trgdists] = make_CF_predictions(...
    session_dir,subject_name,runNum,hemi,trgROI,DoG,seedSig1,seedSig2,seedSig3,seedSig4,TR,SUBJECTS_DIR)

%% Find bold run directories
d = find_bold(session_dir);
%% Set up defaults
if ~exist('hemi','var')
    hemi = 'lh';
end
if ~exist('trgROI','var')
    trgROI='prf_V1';
end
if ~exist('DoG','var');
    DoG = 1;
end
if ~exist('seedSig1','var');
    %seedSig1 = (.5:.5:10)';
    %seedSig1 = (0.5:.5:5)';
    %seedSig1 = (1:1:10)';
    seedSig1 = (0.5:0.5:15)'; % mm cortex
end
if ~exist('seedSig2','var');
    %seedSig2 = (1:.5:3)';
    %seedSig2 = (1:.5:3)';
    %seedSig2 = (1.25:.5:3)';
    seedSig2 = 0;
end
if ~exist('seedSig3','var');
    %seedSig3 = (0:0.125:1)';
    %seedSig3 = (0:0.25:1)';
    seedSig3 = 1;
end
if ~exist('seedSig4','var');
    %seedSig4 = (0:0.125:1)';
    %seedSig4 = (0:0.25:1)';
    seedSig4 = 0;
end
seedSig = {seedSig1 seedSig2 seedSig3 seedSig4};
if ~exist('TR','var')
    TR = 2; % assumes 2s TR
end
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
%% Find bold run directories
d = find_bold(session_dir);
nruns = length(d);
disp(['Session_dir = ' session_dir]);
disp(['Found ' num2str(nruns) ' runs']);
disp(['Working on ' d{runNum}]);
%% Set source and target files
cd(fullfile(session_dir,d{runNum}));
% Get source indices
if strcmp(trgROI,'V1')
    areas = load_nifti(fullfile(session_dir,[hemi '.areas.nii.gz']));
    ecc = load_nifti(fullfile(session_dir,[hemi '.ecc.nii.gz']));
    pol = load_nifti(fullfile(session_dir,[hemi '.pol.nii.gz']));
elseif strcmp(trgROI,'prf_V1')
    areas = load_nifti(fullfile(session_dir,[hemi '.areas_pRF.nii.gz']));
    ecc = load_nifti(fullfile(session_dir,[hemi '.ecc_pRF.nii.gz']));
    pol = load_nifti(fullfile(session_dir,[hemi '.pol_pRF.nii.gz']));
end
% Get target indices, ecc, pol for V1
V1trgind = find(areas.vol==1 | areas.vol ==-1);
V1trgecc = ecc.vol(V1trgind);
V1trgpol = pol.vol(V1trgind);
% Get target indices, ecc, pol for V2
V2trgind = find(areas.vol==2 | areas.vol ==-2);
V2trgecc = ecc.vol(V2trgind);
V2trgpol = pol.vol(V2trgind);
% Get target indices, ecc, pol for V3
V3trgind = find(areas.vol==3 | areas.vol ==-3);
V3trgecc = ecc.vol(V3trgind);
V3trgpol = pol.vol(V3trgind);
% Load target surface file
trgfile = fullfile(session_dir,d{runNum},['sdbrf.tf_surf.' hemi '.nii.gz']);
%% Get distances in target ROI
anatdatadir = fullfile(SUBJECTS_DIR,subject_name);
[verts] = freesurfer_read_surf(fullfile(anatdatadir,'surf',[hemi '.sphere']));
% V1
V1trgdists = nan(length(V1trgind));
progBar = ProgressBar(length(V1trgind),'Calculating V1 distances...');
for i = 1:length(V1trgind);
    vtx = V1trgind(i);
    tmp = spherical_distance(subject_name,vtx,verts,hemi);
    V1trgdists(i,:) = tmp(V1trgind);
    if ~mod(i,100);progBar(i);end
end
% V2
V2trgdists = nan(length(V2trgind));
progBar = ProgressBar(length(V2trgind),'Calculating V2 distances...');
for i = 1:length(V2trgind);
    vtx = V2trgind(i);
    tmp = spherical_distance(subject_name,vtx,verts,hemi);
    V2trgdists(i,:) = tmp(V2trgind);
    if ~mod(i,100);progBar(i);end
end
% V3
V3trgdists = nan(length(V3trgind));
progBar = ProgressBar(length(V3trgind),'Calculating V3 distances...');
for i = 1:length(V3trgind);
    vtx = V3trgind(i);
    tmp = spherical_distance(subject_name,vtx,verts,hemi);
    V3trgdists(i,:) = tmp(V3trgind);
    if ~mod(i,100);progBar(i);end
end
%% Load timecourses
disp('Loading timecourses...');
% load files
trg = load_nifti(trgfile);
trgdims = size(trg.vol);
trgtc = reshape(trg.vol,trgdims(1)*trgdims(2)*trgdims(3),trgdims(4))';
% Pull out relevant timecourses
V1trgtc = trgtc(:,V1trgind);
V2trgtc = trgtc(:,V2trgind);
V3trgtc = trgtc(:,V3trgind);
% Set timecourses with very little variation (var<.1) to flat
V1trgtc = set_to_flat(V1trgtc);
V2trgtc = set_to_flat(V2trgtc);
V3trgtc = set_to_flat(V3trgtc);
disp('done.');
%% Set the tc_shift
% postive values move predicted responses earlier in time, negative values
% move later in time
tc_shift = (-1:0.5:5)/TR;
%% Make sigma list
seedX = 1:length(V1trgdists);
seedSig1 = seedSig{1};
seedSig2 = seedSig{2};
seedSig3 = seedSig{3};
seedSig4 = seedSig{4};
[xList] = ndgrid(seedX,seedSig1,seedSig2,seedSig3,seedSig4);
xList = xList(:);
% Create sigma list for each center location
ct = 0;
for s1 = 1:length(seedSig1);
    for s2 = 1:length(seedSig2);
        for s3 = 1:length(seedSig3);
            for s4 = 1:length(seedSig4);
                ct = ct+1;
                singleSigList(ct,:) = [seedSig1(s1),seedSig1(s1)*seedSig2(s2),seedSig3(s3),seedSig4(s4)];
            end
        end
    end
end
% % Replicate this sigma list for all center locations
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
progBar = ProgressBar(length(tc_shift),'Making predictions...');
for t = 1:length(tc_shift)
    for s = 1:sigTasks;
        % Get V1 center and eccentricity
        V1center = xList(sigvals{s});
        V1x = V1trgdists(:,xList(sigvals{s}));
        V1ecc = V1trgecc(xList(sigvals{s}));
        V1pol = V1trgpol(xList(sigvals{s}));
        [X,Y] = pol2cart(V1pol,V1ecc);
        % Get V1 sigma (mm cortex)
        V1sig = sigList(sigvals{s},:);
        % Convert V1 vals to visual angle
        V1ang = V1sig;
        V1ang(:,1) = V1sig(:,1) ./ cortical_mag(V1ecc,'V1');
        V1ang(:,2) = V1sig(:,2) ./ cortical_mag(V1ecc,'V1');
        % Convert to V2 sigma (mm cortex)
        V2sig = V1sig;
        V2sig(:,1) = V1ang(:,1) .* cortical_mag(V1ecc,'V2');
        V2sig(:,2) = V1ang(:,2) .* cortical_mag(V1ecc,'V2');
        % Convert to V3 sigma (mm cortex)
        V3sig = V1sig;
        V3sig(:,1) = V1ang(:,1) .* cortical_mag(V1ecc,'V3');
        V3sig(:,2) = V1ang(:,2) .* cortical_mag(V1ecc,'V3');
        % Find V2 center
        clear V2center
        [V2allx,V2ally] = pol2cart(V2trgpol,V2trgecc);
        for i = 1:length(X)
            dists = sqrt( (V2allx - X(i)).^2 + (V2ally - Y(i)).^2 );
            [~,V2center(i)] = min(dists);
        end
        V2x = V2trgdists(:,V2center);
        % Find V3 center
        clear V3center
        [V3allx,V3ally] = pol2cart(V2trgpol,V2trgecc);
        for i = 1:length(X)
            dists = sqrt( (V3allx - X(i)).^2 + (V3ally - Y(i)).^2 );
            [~,V3center(i)] = min(dists);
        end
        V3x = V3trgdists(:,V3center);
        % Make predictions for V1
        [~,V1tmpPredResp] = calc_CF_pred(V1x,V1sig,V1trgtc,[],DoG);
        shiftPredResp(1,t,:,sigvals{s}) = shift_tc(V1tmpPredResp,tc_shift(t));
        % Make predictions for V2 (with closest V1 visual field center)
        [~,V2tmpPredResp] = calc_CF_pred(V2x,V2sig,V2trgtc,[],DoG);
        shiftPredResp(2,t,:,sigvals{s}) = shift_tc(V2tmpPredResp,tc_shift(t));
        % Make predictions for V3 (with closest V1 visual field center)
        [~,V3tmpPredResp] = calc_CF_pred(V3x,V3sig,V3trgtc,[],DoG);
        shiftPredResp(3,t,:,sigvals{s}) = shift_tc(V3tmpPredResp,tc_shift(t));
    end
    progBar(t);
end
%% Save file
save(fullfile(session_dir,d{runNum},['CF_predictions.' hemi '.' trgROI '.mat']),'shiftPredResp','tc_shift','xList','sigList','V1trgdists','V2trgdists','V3trgdists','-v7.3');