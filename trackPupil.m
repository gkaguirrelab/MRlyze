function [pupil,glint] = trackPupil(params)

% Tracks the pupil using an input video file, write out an .avi video
%
%   Usage:
%       [pupil,glint]       = trackPupil(params)
%
%   Required input:
%       params.inVideo      = '/path/to/inputFile';
%       params.outVideo     = '/path/to/outVideo.avi';
%       params.outData      = '/path/to/outData.mat'
%
%   Outputs:
%       pupil.X             = X coordinate of pupil center (pixels)
%       pupil.Y             = Y coordinate of pupil center (pixels)
%       pupil.size          = radius of pupil (pixels)
%       glint.X             = X coordinate of glint center (pixels)
%       glint.Y             = Y coordinate of glint center (pixels)
%       glint.size          = radius of glint (pixels)
%
%   Defaults:
%       params.rangeAdjust  = 0.05;         % radius change (+/-) allowed from the previous frame
%       params.threshVals   = [50 175];     % grayscale threshold values for pupil and glint, respectively
%       params.imageSize    = [300 400];    % used to resize input video
%       params.pupilRange   = [10 100];     % initial pupil size range
%       params.glintRange   = [10 15];      % constant glint size range
%       params.glintOut     = 0.1;          % proportion outside of pupil glint is allowed to be. Higher = more outside
%       params.sensitivity  = 0.99;         % [0 1] - sensitivity for 'imfindcircles'. Higher = more circles found
%       params.dilateGlint  = 6;            % used to dialate glint. Higher = more dilation.
%
%   Written by Andrew S Bock Sep 2016

%% set defaults
if ~isfield(params,'rangeAdjust');
    params.rangeAdjust  = 0.05;
end
if ~isfield(params,'binThresh');
    params.threshVals   = [50 175]; % bin for pupil and glint, respectively
end
if ~isfield(params,'imageSize');
    params.imageSize    = [300 400];
end
if ~isfield(params,'pupilRange');
    params.pupilRange   = [10 100];
end
if ~isfield(params,'glintRange');
    params.glintRange   = [10 15];
end
if ~isfield(params,'glintOut');
    params.glintOut     = 0.1;
end
if ~isfield(params,'sensitivity');
    params.sensitivity  = 0.99;
end
if ~isfield(params,'dilateGlint');
    params.dilateGlint  = 6;
end
% Create filter parameters
filtSize                = round([0.01*min(params.imageSize) 0.01*min(params.imageSize) 0.01*min(params.imageSize)]);
% Useful:
%   figure;imshow(I);
%   d = imdistline;
%   % check size of pupil or glint
%   delete(d);
%% Load video
disp('Loading video file, may take a couple minutes...');
inObj                   = VideoReader(params.inVideo);
numFrames               = floor(inObj.Duration*inObj.FrameRate);
grayI                   = zeros([params.imageSize,numFrames],'uint8');
% Convert to gray, resize
for i = 1:numFrames
    thisFrame           = readFrame(inObj);
    tmp                 = rgb2gray(thisFrame);
    grayI(:,:,i)        = imresize(tmp,params.imageSize);
end
clear RGB inObj
outObj                  = VideoWriter(params.outVideo);
open(outObj);
%% Initialize pupil and glint structures
pupilRange              = params.pupilRange;
glintRange              = params.glintRange;
pupil.X                 = nan(numFrames,1);
pupil.Y                 = nan(numFrames,1);
pupil.size              = nan(numFrames,1);
glint.X                 = nan(numFrames,1);
glint.Y                 = nan(numFrames,1);
glint.size              = nan(numFrames,1);
% structuring element to dialate the glint
se                      = strel('disk',params.dilateGlint);
%% Track
progBar = ProgressBar(numFrames,'makingMovie');
ih = fullFigure;
for i = 1:numFrames
    % Get the frame
    I                   = squeeze(grayI(:,:,i));
    % Show the frame
    imshow(I);
    % Filter for pupil
    padP                = padarray(I,[size(I,1)/2 size(I,2)/2], 128);
    h                   = fspecial('gaussian',[filtSize(1) filtSize(2)],filtSize(3));
    pI                  = imfilter(padP,h);
    pI = pI(size(I,1)/2+1:size(I,1)/2+size(I,1),size(I,2)/2+1:size(I,2)/2+size(I,2));
    % Binarize pupil
    binP                = ones(size(pI));
    binP(pI<params.threshVals(1))   = 0;
    % Filter for glint
    gI                  = ones(size(I));
    gI(I<params.threshVals(2)) = 0;
    padG                = padarray(gI,[size(I,1)/2 size(I,2)/2], 0);
    h                   = fspecial('gaussian',[filtSize(1) filtSize(2)],filtSize(3));
    gI                  = imfilter(padG,h);
    gI = gI(size(I,1)/2+1:size(I,1)/2+size(I,1),size(I,2)/2+1:size(I,2)/2+size(I,2));
    % Binarize glint
    binG                = zeros(size(gI));
    binG(gI>0.01)       = 1;
    dbinG               = imdilate(binG,se);
    % Find the pupil
    [pCenters, pRadii]  = imfindcircles(binP,pupilRange,'ObjectPolarity','dark',...
        'Sensitivity',params.sensitivity);
    % Find the glint
    [gCenters, gRadii]  = imfindcircles(dbinG,glintRange,'ObjectPolarity','bright',...
        'Sensitivity',params.sensitivity);
    % Remove glints outside the pupil
    if ~isempty(pCenters) && ~isempty(gCenters)
        dists           = sqrt( (gCenters(:,1) - pCenters(1,1)).^2 + (gCenters(:,2) - pCenters(1,2)).^2 );
        gCenters(dists>(1 + params.glintOut)*(pRadii(1)),:) = [];
        gRadii(dists>(1 + params.glintOut)*(pRadii(1))) = [];
    end
    % Visualize the pupil and glint on the image
    if ~isempty(pCenters) && ~isempty(gCenters)
        pupil.X(i)      = pCenters(1,1);
        pupil.Y(i)      = pCenters(1,2);
        pupil.size(i)   = pRadii(1);
        glint.X(i)      = gCenters(1,1);
        glint.Y(i)      = gCenters(1,2);
        glint.size(i)   = gRadii(1);
        viscircles(pCenters(1,:),pRadii(1),'Color','r');
        viscircles(gCenters(1,:),gRadii(1),'Color','b');
        pupilRange(1)   = min(floor(pRadii(1)*(1-params.rangeAdjust)),params.pupilRange(2));
        pupilRange(2)   = max(ceil(pRadii(1)*(1 + params.rangeAdjust)),params.pupilRange(1));
    end
    % Adjust range if pupil is not found
    if ~isempty(pCenters)
        pupilRange(1)   = max(ceil(pupilRange(1)*(1 - params.rangeAdjust)),params.pupilRange(1));
        pupilRange(2)   = min(ceil(pupilRange(2)*(1 + params.rangeAdjust)),params.pupilRange(2));
    end
    frame               = getframe(ih);
    writeVideo(outObj,frame);
    if ~mod(i,10);progBar(i);end;
end
close(outObj);
close(ih);
save(params.outData,'pupil','glint');