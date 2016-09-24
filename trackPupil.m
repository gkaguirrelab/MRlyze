function trackPupil(params)

% Tracks the pupil using an input video file, write out an .avi video
%
%   Usage:
%       trackPupil(params)
%
%   Required input:
%       params.inVideo      = '/path/to/inputFile';
%       params.outVideo     = '/path/to/outputFile';
%
%   Defaults:
%       params.rangeAdjust  = 0.15;         % radius change (+/-) allowed from the previous frame
%       params.threshVals   = [50 175];     % grayscale threshold values for pupil and glint, respectively
%       params.imageSize    = [900 1200];   % used to resize input video
%       params.pupilRange   = [50 300];     % initial pupil size range
%       params.glintRange   = [10 30];      % constant glint size range
%
%   Written by Andrew S Bock Sep 2016

%% set defaults
if ~isfield(params,'rangeAdjust');
    params.rangeAdjust  = 0.15;
end
if ~isfield(params,'binThresh');
    params.threshVals   = [50 175]; % bin for pupil and glint, respectively
end
if ~isfield(params,'imageSize');
    params.imageSize    = [900 1200];
end
if ~isfield(params,'pupilRange');
    params.pupilRange   = [50 300];
end
if ~isfield(params,'glintRange');
    params.glintRange   = [10 30];
end
% Create filter parameters
filtSize                = [0.05*min(params.imageSize) 0.05*min(params.imageSize) 0.01*min(params.imageSize)];

% Useful:
%   figure;imshow(I);
%   d = imdistline;
%   % check size of pupil or glint
%   delete(d);
%% Load video
disp('Loading video file, may take a couple minutes...');
inObj                   = VideoReader(params.inVideo);
%NumberOfFrames         = inObj.NumberOfFrames;
%RGB                     = read(inObj);
NumberOfFrames          = floor(inObj.Duration*inObj.FrameRate);
grayI                   = zeros([params.imageSize,NumberOfFrames],'uint8');
% Convert to gray, resize
for i = 1:NumberOfFrames
    thisFrame           = readFrame(inObj);
    tmp                 = rgb2gray(thisFrame);
    grayI(:,:,i)        = imresize(tmp,params.imageSize);
end
clear RGB inObj
outObj                  = VideoWriter(params.outVideo);
open(outObj);
%% Track
progBar = ProgressBar(NumberOfFrames,'makingMovie');
ih = fullFigure;
for i = 1:NumberOfFrames
    % Get the frame
    I                   = squeeze(grayI(:,:,i));
    % Filter for glint
    gI                  = ones(size(I));
    gI(I<params.threshVals(2)) = 0;
    padG                = padarray(gI,[size(I,1)/2 size(I,2)/2], 0);
    h                   = fspecial('gaussian',[filtSize(1) filtSize(2)],filtSize(3));
    gI                  = imfilter(padG,h);
    gI = gI(size(I,1)/2+1:size(I,1)/2+size(I,1),size(I,2)/2+1:size(I,2)/2+size(I,2));
    % Filter for pupil
    padP                = padarray(I,[size(I,1)/2 size(I,2)/2], 128);
    h                   = fspecial('gaussian',[filtSize(1) filtSize(2)],filtSize(3));
    fI                  = imfilter(padP,h);
    fI = fI(size(I,1)/2+1:size(I,1)/2+size(I,1),size(I,2)/2+1:size(I,2)/2+size(I,2));
    % Binarize the image
    binP                = ones(size(fI));
    binP(fI<params.threshVals(1))   = 0;
    binG                = zeros(size(gI));
    binG(gI>0.1)        = 1;
    % Show the frame
    imshow(I);
    % Find the pupil
    [pCenters, pRadii]  = imfindcircles(binP,params.pupilRange,'ObjectPolarity','dark','Sensitivity',0.99);
    % Find the glint
    [gCenters, gRadii]  = imfindcircles(binG,params.glintRange,'ObjectPolarity','bright','Sensitivity',0.99);
    % Remove glints outside the pupil
    if ~isempty(pCenters) && ~isempty(gCenters)
        dists = sqrt( (gCenters(:,1) - pCenters(1,1)).^2 + (gCenters(:,2) - pCenters(1,2)).^2 );
        gCenters(dists>pRadii(1),:) = [];
        gRadii(dists>pRadii(1)) = [];
    end
    % Visualize the pupil and glint on the image
    if ~isempty(pCenters) && ~isempty(gCenters)
        viscircles(pCenters(1,:),pRadii(1),'Color','r');
        viscircles(gCenters(1,:),gRadii(1),'Color','b');
        params.pupilRange = [floor(pRadii(1)*(1-params.rangeAdjust)) ceil(pRadii(1)*(1 + params.rangeAdjust))];
    end
    frame               = getframe(ih);
    writeVideo(outObj,frame);
    progBar(i);
end
close(outObj);