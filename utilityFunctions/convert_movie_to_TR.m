function convert_movie_to_TR(saveFile,movieFile,movieBlock,TR)

% Convert RGB movie with a certain frames per second to grayscale, sampled
%   at the TR
%
%   Written by Andrew S Bock Jul 2015


%% Set defaults
if ~exist('saveFile','var')
   error('No ''saveFile'' input provided');
end
if ~exist('movieFile','var')
    movieFile = '/jet/abock/data/Retinotopy/ALL.mov';
end
if ~exist('movieBlock','var')
    movieBlock = [0 600]; % [start stop] of movie, in seconds
end
if ~exist('TR','var')
    TR = 2; % repetition time (seconds)
end
%%   Convert movie to grayscale and average images by TR
vid = VideoReader(movieFile);
frameDur = movieBlock*vid.FrameRate; % [start stop] of movie, in frames
numFrames = frameDur(2) - frameDur(1);
numTRs = (movieBlock(2) - movieBlock(1))/TR;
[framesPerTR] = calc_tasks(numFrames,numTRs);
mov = zeros(vid.Height,vid.Width,numTRs,'uint8');
ct = 0;
progBar = ProgressBar(numTRs,'Looping through movie...');
for tt = 1:numTRs
    tmpMOV = zeros(vid.Height,vid.Width,framesPerTR(tt),'uint8');
    for kk = 1:framesPerTR(tt);
        ct = ct + 1;
        tmpRGB = read(vid,frameDur(1) + ct);
        tmpMOV(:,:,kk) = rgb2gray(tmpRGB);
    end
    mov(:,:,tt) = mean(tmpMOV,3);
    progBar(tt);
end
%% Save converted movie file
save(saveFile,'mov');