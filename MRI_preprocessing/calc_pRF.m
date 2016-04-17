function [prfs] = calc_pRF(srctc,TR,nFrames,fieldSize,imFile,paramsFile)

%% Set hard coded defaults
% Repetition time
disp(['TR = ' num2str(TR)]);
% Time to peak
peakt = 3:.5:10; % time to HRF peak 
disp(['peak HRF values = ' num2str(peakt)]);
% number of TRs in scan
%nFrames = 300;
disp(['nFrames = ' num2str(nFrames)]);
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
% visual angle radius (UPenn 7T = 6.2116)
%fieldSize = 6.2116 ;
% allowable pRF centers
gridscale = 2; % double the size of the visual field based on the stimulated region
maxXY = gridscale*fieldSize;
% based on mrVista (default is 50)
GridPoints = 50;
% sample rate in visual angle
sampleRate = maxXY./GridPoints;
%% Load image and parameter files
disp('Loading images and parameters...');
I = load(imFile);
P = load(paramsFile);
% Binarize images (i.e. 1 if image, 0 if background)
bk = P.params.display.backColorIndex; % find 'black' value
I.images = I.images ~= bk;
disp('done.');
%% Convert srctc to single
srctc = single(srctc);

%% Make sigma list
tmps = ctrSigbnd(1):ctrSigbnd(2):Sigborder;
linsigs = tmps(1:end-1);
logsigs = logspace(log10(Sigborder),log10(ctrSigbnd(3)),nLogSigs);
ctrSigma = [linsigs logsigs];
%% Create full sigma list
% Create sigma list for each center location
clear sigList
ct = 0;
for s1 = 1:length(ctrSigma);
    for s2 = 1:length(srdSigma);
        for s3 = 1:length(ctrBetas);
            for s4 = 1:length(srdBetas);
                ct = ct+1;
                sigList(ct,:) = [ctrSigma(s1),ctrSigma(s1)*srdSigma(s2),ctrBetas(s3),srdBetas(s4)];
            end
        end
    end
end
badind = sigList(:,3) == 0 & sigList(:,4) == 0;
sigList(badind,:) = [];
%% Create grid of visual space (x,y) with corresponding sigma values (z)
mygrid = -maxXY:sampleRate:maxXY;
[tmpx0,tmpy0]=meshgrid(mygrid,mygrid);
tmpx0 = tmpx0(:);
tmpy0 = tmpy0(:);
% Save visual field x, y, and sigma values
% define parameters
stim.x0 = repmat(tmpx0,size(sigList,1),1);
stim.y0 = repmat(tmpy0,size(sigList,1),1);
stim.sigs = repmat(sigList,size(tmpx0,1),1);
%% Create sampling grid for images
disp('Creating sampling grid for images...');
mygrid = -fieldSize:sampleRate:fieldSize;
[X,Y]=meshgrid(mygrid,mygrid);
X = X(:);
Y = Y(:);
% Downsample images to sampling grid
nImages = size(I.images, 3);
resampled = zeros(length(X),nImages);
for ii = 1:nImages
    tmp_im = imresize(I.images(:,:,ii), 1+2*[GridPoints/gridscale GridPoints/gridscale]);
    resampled(:, ii) = tmp_im(:);
end
% Temporally downsample to 1 image per TR (by averaging filtered images)
% length of 1 TR
seq         = P.stimulus.seq;           %index to image number
seqTiming   = P.stimulus.seqtiming;     %time in s for each image onset
% Temporally downsample
images = zeros(length(X),nFrames);
% Specificy the onset time and offset time of each image frame
imOnset = seqTiming;
imOffset = shift(seqTiming, -1);
imOffset(end) = TR * nFrames;
% Create binary matrix corresponding to where stimulus was present
for f = 1:nFrames
    % specify the onset time and offset time of this TR
    frameOnset  = TR * (f-1);
    frameOffset = TR * f;
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
    img = img / TR;
    images(:,f) = img;
end
images = single(images);
images_bak = images;
% now scale amplitude according to the sample rate, for %BOLD/degree2
images = images.*(sampleRate.^2);
disp('done.');
%% Add black around stimulus region, to model the actual visual field (not just the bars)
mygrid = -maxXY:sampleRate:maxXY;
[x,y]=meshgrid(mygrid,mygrid);
images = reshape(images,sqrt(size(images,1)),sqrt(size(images,1)),size(images,2));
tmpimages = zeros([size(x) size(images,3)]);
diffsize = size(tmpimages,1) - size(images,1);
inds = diffsize/2+1:diffsize/2+size(images,1);
tmpimages(inds,inds,:) = images;
stim.X = x(:);
stim.Y = y(:);
images = reshape(tmpimages,size(tmpimages,1)*size(tmpimages,2),size(tmpimages,3));
%% Make predicted timecoures from stimulus images
nn = numel(stim.x0);
p = [(1:ceil(nn./500):nn-2) nn+1];
maxs = p(numel(p))-1;
prediction = zeros(size(images,2),nn*length(peakt),'single');
% Create HRF
progBar = ProgressBar(length(peakt),'Making predictions...');
for h = 1:length(peakt);
    disp(['Predictions for HRF ' num2str(peakt(h)) 's peak - ' num2str(h) ' of ' num2str(length(peakt))]);
    Hrf = doubleGammaHrf(TR,[peakt(h) peakt(h)+10]);
    %progBar = ProgressBar(numel(s)-1,'Making predictions...');
    for n=1:numel(p)-1
        X = stim.X;
        Y = stim.Y;
        sigs = stim.sigs(p(n):p(n+1)-1,:);
        x0 = stim.x0(p(n):p(n+1)-1);
        y0 = stim.y0(p(n):p(n+1)-1);
        % Allow sigma, x,y to be a matrix so that the final output will be
        % size(X,1) by size(x0,2). This way we can make many RFs at the same time.
        % Here I assume that all parameters are given.
        if numel(sigs)~=1,
            sz1 = size(X,1);
            sz2 = size(sigs,1);
            X   = repmat(X,1,sz2);
            Y   = repmat(Y(:),1,sz2);
            x0 = repmat(x0',sz1,1);
            y0 = repmat(y0',sz1,1);
            sigs = repmat(sigs,1,1,sz1);
            sigs = permute(sigs,[3 1 2]);
        end;
        % Translate grid so that center is at RF center
        nX = X - x0;   % positive x0 moves center right
        nY = Y - y0;   % positive y0 moves center up
        % make gaussian on current grid
        %cG =   sigs(:,:,3).*(exp (-(nY.^2 + nX.^2) ./ (2*sigs(:,:,1).^2)));
        %sG =   sigs(:,:,4).*(exp (-(nY.^2 + nX.^2) ./ (2*sigs(:,:,2).^2)));
        %rf = cG - sG;
        rf = exp (-(nY.^2 + nX.^2) ./ (2*sigs(:,:,1).^2));
        % Convolve images with Hrf
        stim.allstimimages = filter(Hrf,1, images');
        % Convolve images (with Hrf) with receptive field
        pred = stim.allstimimages*rf;
        % Set timecourses with very little variation (var<0.1) to flat
        % above images are set to be %BOLD/degree2
        pred = set_to_flat(pred);
        % store
        prediction(:,(p(n)+(h-1)*maxs:p(n+1)-1+(h-1)*maxs)) = pred;
        %if ~mod(n,round((numel(s)-1)/20));progBar(n);end
        %progBar(n);
    end
    progBar(h);
end
% Create x0, y0, sig, and peak matrices (for multiple peak times)
stimx0mat = repmat(stim.x0,[length(peakt) 1]);
stimy0mat = repmat(stim.y0,[length(peakt) 1]);
stimsigmat = repmat(stim.sigs,[length(peakt) 1]);
peaktmat = repmat(peakt,[length(stimx0mat)/length(peakt) 1]);
%% Remove the predictions with very little variation (i.e. flat timecourses)
predvar = zeros(1,size(prediction,2));
for i = 1:size(prediction,2);
    predvar(i) = var(prediction(:,i));
end
badind = predvar<0.1;
prediction(:,badind) = [];
stimx0mat(badind) = [];
stimy0mat(badind) = [];
stimsigmat(badind,:) = [];
peaktmat(badind) = [];
%% Break up sigma list matrix
[predPerTask,predTasks] = calc_tasks(size(prediction,2),ceil(size(prediction,2)/1000));
predidx = [];
for i = 1:predTasks
    if isempty(predidx);
        predidx = [1,predPerTask(i)];
    else
        predidx = [predidx;[predidx(end,2)+1,predidx(end,2)+predPerTask(i)]];
    end
    predvals{i} = predidx(i,1):predidx(i,2);
end
%% Initialize variables
prfs.co = zeros(size(srctc,2),1,'single');
prfs.cox = zeros(size(srctc,2),1,'single');
prfs.coy = zeros(size(srctc,2),1,'single');
prfs.cosig = zeros(size(srctc,2),4,'single');
prfs.copeakt = zeros(size(srctc,2),1,'single');
prfs.var = zeros(size(srctc,2),1,'single');
prfs.var_as_co = zeros(size(srctc,2),1,'single');
prfs.varx = zeros(size(srctc,2),1,'single');
prfs.vary = zeros(size(srctc,2),1,'single');
prfs.varsig = zeros(size(srctc,2),4,'single');
prfs.varpeakt = zeros(size(srctc,2),1,'single');
%% Create Predicted Responses
%disp('Finding Connective Fields...');
%poolobj = gcp; % start parpool
disp('Finding Population Receptive Fields...');
tstart = clock;
ProgressBar_parfor(tstart);
%progBar = ProgressBar(sigTasks,'Finding Connective Fields...');
for p = 1:predTasks;
    tmp_co = corr(srctc,prediction(:,predvals{p}));
    tmp_co = single(tmp_co);
    % Save best co
    [prfsco,prfscoseed] = max(tmp_co,[],2);
    newcoind = prfs.co<prfsco;
    prfs.co(newcoind) = prfsco(newcoind);
    prfs.cox(newcoind) = stimx0mat(predvals{p}(prfscoseed(newcoind)));
    prfs.coy(newcoind) = stimy0mat(predvals{p}(prfscoseed(newcoind)));
    prfs.cosig(newcoind,:) = stimsigmat(predvals{p}(prfscoseed(newcoind)),:);
    prfs.copeakt(newcoind) = peaktmat(predvals{p}(prfscoseed(newcoind)));
    % Save best var
    [prfsvar,prfsvarseed] = max(tmp_co.^2,[],2);
    newvarind = prfs.var<prfsvar;
    prfs.var(newvarind) = prfsvar(newvarind);
    prfs.varx(newvarind) = stimx0mat(predvals{p}(prfsvarseed(newvarind)));
    prfs.vary(newvarind) = stimy0mat(predvals{p}(prfsvarseed(newvarind)));
    prfs.varsig(newvarind,:) = stimsigmat(predvals{p}(prfsvarseed(newvarind)),:);
    prfs.varpeakt(newvarind) = peaktmat(predvals{p}(prfsvarseed(newvarind)));
    vcind = find(newvarind);
    for vc = 1:length(vcind);
        prfs.var_as_co(vcind(vc)) = tmp_co(vcind(vc),prfsvarseed(vcind(vc)));
    end
    ProgressBar_parfor(tstart,p,predTasks);
end
ProgressBar_parfor(tstart,'clean'); % Clean up files and display total loop time
disp('done.');
%% Convert from cartesian to polar
[prfs.copol,prfs.coecc] = cart2pol(prfs.cox,prfs.coy);
[prfs.varpol,prfs.varecc] = cart2pol(prfs.varx,prfs.vary);