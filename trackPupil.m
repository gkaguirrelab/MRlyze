function trackPupil(params)

% Tracks the pupil using an input video file, write out an .avi video
%
%   Usage: 
%       trackPupil(params)
%
%   Written by Andrew S Bock Sep 2016

%% set defaults
if ~isfield(params,'video')
   params.video = '/Users/abock/Dropbox-Aguirre-Brainard-Lab/TOME_data/session1_restandstructure/tome_3001/081916/EyeTracking/rfMRI_REST_AP_run01_raw.mov';
end
if ~isfield(params,'NumberOfFrames');
    params.NumberOfFrames = 100;
end
rangeAdjust = 0.1;
imageSize = [900 1200];
filtSize = [40 40 10];
glintRange = [10 20];
pupilRange = [100 200];
% Good videos:
% '/Users/abock/Dropbox-Aguirre-Brainard-Lab/TOME_data/session1_restAndStructure/TOME_3004/091916/EyeTracking/rfMRI_REST_AP_run01_raw.mov';
% Bad videos:
% '/Users/abock/Dropbox-Aguirre-Brainard-Lab/TOME_data/session1_restandstructure/tome_3001/081916/EyeTracking/rfMRI_REST_AP_run01_raw.mov';
% '/Users/abock/Dropbox-Aguirre-Brainard-Lab/TOME_data/session2_spatialstimuli/tome_3002/082616/EyeTracking/tfMRI_MOVIE_AP_run02_raw.mov';
%% Load video
obj = VideoReader(params.video);
NumberOfFrames = obj.NumberOfFrames;
aviobj = VideoWriter('BockPupilTracking');
open(aviobj);
%% Track
progBar = ProgressBar(NumberOfFrames,'makingMovie');
ih = fullFigure;
for i = 1:NumberOfFrames
    RGB                 = read(obj,i);
    I                   = rgb2gray(RGB);
    I                   = imresize(I,imageSize);
    h                   = fspecial('gaussian',[filtSize(1) filtSize(2)],filtSize(3));
    fI                  = imfilter(I,h); % smooth image
    % Pupil (smoothed image)
    [pCenters, pRadii]  = imfindcircles(fI,pupilRange,'ObjectPolarity','dark','Sensitivity',0.99);
    % Glint (unsmoothed image)
    [gCenters, gRadii]  = imfindcircles(fI,glintRange,'ObjectPolarity','bright','Sensitivity',0.99);
    imshow(I);
    if ~isempty(pCenters) && ~isempty(gCenters)
        ph              = viscircles(pCenters(1,:),pRadii(1),'Color','r');
        gh              = viscircles(gCenters(1,:),gRadii(1),'Color','b');
        pupilRange      = [floor(pRadii(1)*(1-rangeAdjust)) ceil(pRadii(1)*(1 + rangeAdjust))];
        glintRange      = [floor(gRadii(1)*(1-rangeAdjust)) ceil(gRadii(1)*(1 + rangeAdjust))];
    end
    frame               = getframe(ih);
    writeVideo(aviobj,frame);
    progBar(i);
end
close(aviobj);