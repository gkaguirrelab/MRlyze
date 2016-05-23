function createPredpRF(imFile,paramsFile,outFile,params)

% Creates predicted pRF time-series data
%
%   Usage:
%   createPredpRF(imFile,paramsFile,outFile,params)
%
%   Defaults:
%     params.fieldSize         = 6.2116; % radius of stimuluated visual field in degrees visual angle
%     params.TR                = 2; % TR is seconds
%     params.peakHRF           = 3:.5:10; % time to HRF peak
%     params.nFrames           = 300; % number of images
%     params.cntrSig           = [1/4 1/4 20]; % min/lin_step/max sigma in degrees visual angle
%     params.sigBorder         = 5; % Linear scale for visual field less than 5 degrees
%     params.nLogSigs          = 6; % number of log sigmas
%     params.srndSig           = 0; % scale of surround sigma 
%     params.ctrBetas          = 1; % center amplitude 
%     params.srndBetas         = 0; % surround amplitude 
%     params.gridScale         = 2; % scale the size of the stimulated visual field
%
%   Outputs:
%   'predTCs','stimX0','stimy0','stimsig','peakHRF'
%
%   Written by Andrew S Bock May 2016

%% Set defaults
if ~exist('params','var')
    params.fieldSize         = 6.2116; % radius of stimuluated visual field in degrees visual angle
    params.TR                = 2; % TR is seconds
    params.peakHRF           = 3:.5:10; % time to HRF peak
    params.nFrames           = 300; % number of images
    params.cntrSig           = [1/4 1/4 20]; % min/lin_step/max sigma in degrees visual angle
    params.sigBorder         = 5; % Linear scale for visual field less than 5 degrees
    params.nLogSigs          = 6; % number of log sigmas
    params.srndSig           = 0; % scale of surround sigma % (1.5:.75:3)';
    params.ctrBetas          = 1; % center amplitude % (0:.25:1)';
    params.srndBetas         = 0; % surround amplitude % (0:.25:1)';
    params.gridScale         = 2; % scaling of the stimulated visual field
end
maxXY                       = params.gridScale*params.fieldSize;
GridPoints                  = 50; % based on mrVista (default is 50)
sampleRate                  = maxXY./GridPoints; % sample rate in visual angle
%% Load image and parameter files
I = load(imFile);
P = load(paramsFile);
% Binarize images (i.e. 1 if image, 0 if background)
bk = P.params.display.backColorIndex; % find 'black' value
I.images = I.images ~= bk;
%% Make sigma list
tmps = params.cntrSig(1):params.cntrSig(2):params.sigBorder;
linsigs = tmps(1:end-1);
logsigs = logspace(log10(params.sigBorder),log10(params.cntrSig(3)),params.nLogSigs);
sigmaVals = [linsigs logsigs];
%% Create full sigma list
% Create sigma list for each center location
clear sigList
ct = 0;
for s1 = 1:length(sigmaVals);
    for s2 = 1:length(params.srndSig);
        for s3 = 1:length(params.ctrBetas);
            for s4 = 1:length(params.srndBetas);
                ct = ct+1;
                sigList(ct,:) = [sigmaVals(s1),sigmaVals(s1)*params.srndSig(s2),...
                    params.ctrBetas(s3),params.srndBetas(s4)];
            end
        end
    end
end
badind = sigList(:,3) == 0 & sigList(:,4) == 0;
sigList(badind,:) = [];
%% Create grid of visual space (x,y) with corresponding sigma values (z)
tmpgrid = -maxXY:sampleRate:maxXY;
[tmpx0,tmpy0]=meshgrid(tmpgrid,tmpgrid);
tmpx0 = tmpx0(:);
tmpy0 = tmpy0(:);
% Save visual field x, y, and sigma values
% define parameters
stim.x0 = repmat(tmpx0,size(sigList,1),1);
stim.y0 = repmat(tmpy0,size(sigList,1),1);
stim.sigs = repmat(sigList,size(tmpx0,1),1);
%% Create sampling grid for images
tmpgrid = -params.fieldSize:sampleRate:params.fieldSize;
[X,~]=meshgrid(tmpgrid,tmpgrid);
X = X(:);
%% Downsample images to sampling grid
nImages = size(I.images, 3);
resampled = zeros(length(X),nImages);
for ii = 1:nImages
    tmp_im = imresize(I.images(:,:,ii), 1+2*[GridPoints/params.gridScale GridPoints/params.gridScale]);
    resampled(:, ii) = tmp_im(:);
end
%% Temporally downsample to 1 image per TR (by averaging filtered images)
% length of 1 TR
seq         = P.stimulus.seq;           %index to image number
seqTiming   = P.stimulus.seqtiming;     %time in s for each image onset
% Temporally downsample
images = zeros(length(X),params.nFrames);
%% Specificy the onset time and offset time of each image frame
imOnset = seqTiming;
imOffset = shift(seqTiming, -1);
imOffset(end) = params.TR * params.nFrames;
%% Create binary matrix corresponding to where stimulus was present
for f = 1:params.nFrames
    % specify the onset time and offset time of this TR
    frameOnset  = params.TR * (f-1);
    frameOffset = params.TR * f;
    % calculate the temporal overlap of each image frame with this TR
    imDur = min(imOffset, frameOffset) - max(imOnset, frameOnset);
    imDur = max(imDur, 0);
    imDur(imDur < .001) = 0;
    img = zeros(size(resampled(:,1)));
    % weighted average (only loop over seq we need)
    ii = find(imDur>0); ii = ii(:)';
    for im = ii
        img = img + imDur(im) * resampled(:, seq(im));
    end
    img = img / params.TR;
    images(:,f) = img;
end
images = single(images);
% now scale amplitude according to the sample rate, for %BOLD/degree2
images = images.*(sampleRate.^2);
%% Add black around stimulus region, to model the actual visual field (not just the bars)
tmpgrid = -maxXY:sampleRate:maxXY;
[x,y]=meshgrid(tmpgrid,tmpgrid);
images = reshape(images,sqrt(size(images,1)),sqrt(size(images,1)),size(images,2));
tmpimages = zeros([size(x) size(images,3)]);
diffsize = size(tmpimages,1) - size(images,1);
inds = diffsize/2+1:diffsize/2+size(images,1);
tmpimages(inds,inds,:) = images;
stim.X = x(:);
stim.Y = y(:);
images = reshape(tmpimages,size(tmpimages,1)*size(tmpimages,2),size(tmpimages,3));
%% Break up into smaller matrices
nn = numel(stim.x0); % grid points
[predPerTask,predTasks] = calc_tasks(nn,ceil(nn/1000));
predidx = [];
for i = 1:predTasks
    if isempty(predidx);
        predidx = [1,predPerTask(i)];
    else
        predidx = [predidx;[predidx(end,2)+1,predidx(end,2)+predPerTask(i)]];
    end
    predvals{i} = predidx(i,1):predidx(i,2);
end
%% Make predicted timecoures from stimulus images
predTCs = zeros(size(images,2),nn*length(params.peakHRF),'single');
% Create HRF
progBar = ProgressBar(length(params.peakHRF),'Making predictions...');
for h = 1:length(params.peakHRF);
    disp(['Predictions for HRF ' num2str(params.peakHRF(h)) 's peak - ' num2str(h) ' of ' num2str(length(params.peakHRF))]);
    HRF = doubleGammaHrf(params.TR,[params.peakHRF(h) params.peakHRF(h)+10]);
    for n=1:length(predvals)
        X = stim.X;
        Y = stim.Y;
        sigs = stim.sigs(predvals{n},:);
        x0 = stim.x0(predvals{n});
        y0 = stim.y0(predvals{n});
        % Allow x, y, sigma to be a matrix so that the final output will be
        % size(X,1) by size(x0,2). This way we can make many RFs at the same time.
        if numel(sigs)~=1,
            sz1 = size(stim.X,1);
            sz2 = size(sigs,1);
            X   = repmat(X,1,sz2);
            Y   = repmat(Y(:),1,sz2);
            x0 = repmat(x0',sz1,1);
            y0 = repmat(y0',sz1,1);
            sigs = repmat(sigs,1,1,sz1);
            sigs = permute(sigs,[3 1 2]);
        end
        % Translate grid so that center is at RF center
        nX = X - x0;   % positive x0 moves center right
        nY = Y - y0;   % positive y0 moves center up
        % make gaussian on current grid
        rf = exp (-(nY.^2 + nX.^2) ./ (2*sigs(:,:,1).^2));
        % Convolve images with HRF
        stim.allstimimages = filter(HRF,1, images');
        % Convolve images (with HRF) with receptive field
        pred = stim.allstimimages*rf;
        % Set timecourses with very little variation (var<0.1) to flat
        % above images are set to be %BOLD/degree2
        pred = set_to_flat(pred);
        % store
        predTCs(:,(predvals{n} + (h-1)*nn)) = pred;
    end
    progBar(h);
end
% Create x0, y0, sig, and peak matrices (for multiple peak times)
stimX0 = repmat(stim.x0,[length(params.peakHRF) 1]);
stimY0 = repmat(stim.y0,[length(params.peakHRF) 1]);
stimSig = repmat(stim.sigs,[length(params.peakHRF) 1]);
peakHRF = repmat(params.peakHRF,[length(stimX0)/length(params.peakHRF) 1]);
%% Remove the predictions with very little variation (i.e. flat timecourses)
predvar = zeros(1,size(predTCs,2));
for i = 1:size(predTCs,2);
    predvar(i) = var(predTCs(:,i));
end
badind = predvar<0.1;
predTCs(:,badind) = [];
stimX0(badind) = [];
stimY0(badind) = [];
stimSig(badind,:) = [];
peakHRF(badind) = [];
%% Save output
save(outFile,'predTCs','stimX0','stimY0','stimSig','peakHRF','v-7.3');