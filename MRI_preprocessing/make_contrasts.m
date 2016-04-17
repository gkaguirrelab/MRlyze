function [contrasts] = make_contrasts(runNum,outputDir,blockdur,TR,maketext)

% Creates a Nx3 matrix, with N number of runs, and 3 conditions per run.
% The duration of each block (sec) is defined by 'blockdur'.
%
%   Usage:
%   make_contrasts(runs,blockdur,maketext,outname)
%
%   defaults:
%   runNum = 1;
%   blockdur = 16 % in seconds
%   TR = 2 % in seconds
%   maketext = 1% print text files, e.g. for use with FSL
%   outputDir = no default, must specify the output name of the text file
%
% The occurance of each condition is based on the De Bruijn sequence for
% that run:
%
% DeBruijn_sequence(1).vals = [1,2,1,1,3,3,2,2,3,1,2,1,1,3,3,2,2,3,1];
% DeBruijn_sequence(2).vals = [1,3,3,2,1,1,2,2,3,1,3,3,2,1,1,2,2,3,1];
% DeBruijn_sequence(3).vals = [1,2,2,3,1,3,3,2,1,1,2,2,3,1,3,3,2,1,1];
% DeBruijn_sequence(4).vals = [1,3,2,3,3,1,2,2,1,1,3,2,3,3,1,2,2,1,1];
% DeBruijn_sequence(5).vals = [1,3,2,2,3,3,1,1,2,1,3,2,2,3,3,1,1,2,1];
% DeBruijn_sequence(6).vals = [1,3,1,1,2,2,3,3,2,1,3,1,1,2,2,3,3,2,1];
% DeBruijn_sequence(7).vals = [1,2,2,3,2,1,3,3,1,1,2,2,3,2,1,3,3,1,1];
% DeBruijn_sequence(8).vals = [1,2,3,3,2,2,1,3,1,1,2,3,3,2,2,1,3,1,1];
% DeBruijn_sequence(9).vals = [1,3,3,2,1,2,2,3,1,1,3,3,2,1,2,2,3,1,1];
% DeBruijn_sequence(10).vals = [1,2,2,1,3,2,3,3,1,1,2,2,1,3,2,3,3,1,1];
% DeBruijn_sequence(999).vals = [1,2,1,1,3,3,2,2,3,1];
% if run == 999
%     DeBruijn = DeBruijn_sequence(run).vals;
% elseif run > 10
%     DeBruijn = DeBruijn_sequence(mod(run,10)).vals;
% else
%     DeBruijn = DeBruijn_sequence(run).vals;
% end
%
%   Written by Andrew S Bock Oct 2014
%
% TO DO:
%   - Make this more generic

%% Set defaults
if ~exist('runNum','var')
    runNum = 1;
end
if ~exist('blockdur','var')
    blockdur = 16; % seconds
end
if ~exist('TR','var')
    TR=2; % seconds
end
if ~exist('maketext','var')
    maketext = 1;
end
%% Define De Bruijn Sequences
DeBruijn_sequence(1).vals = [1,2,1,1,3,3,2,2,3,1,2,1,1,3,3,2,2,3,1];
DeBruijn_sequence(2).vals = [1,3,3,2,1,1,2,2,3,1,3,3,2,1,1,2,2,3,1];
DeBruijn_sequence(3).vals = [1,2,2,3,1,3,3,2,1,1,2,2,3,1,3,3,2,1,1];
DeBruijn_sequence(4).vals = [1,3,2,3,3,1,2,2,1,1,3,2,3,3,1,2,2,1,1];
DeBruijn_sequence(5).vals = [1,3,2,2,3,3,1,1,2,1,3,2,2,3,3,1,1,2,1];
DeBruijn_sequence(6).vals = [1,3,1,1,2,2,3,3,2,1,3,1,1,2,2,3,3,2,1];
DeBruijn_sequence(7).vals = [1,2,2,3,2,1,3,3,1,1,2,2,3,2,1,3,3,1,1];
DeBruijn_sequence(8).vals = [1,2,3,3,2,2,1,3,1,1,2,3,3,2,2,1,3,1,1];
DeBruijn_sequence(9).vals = [1,3,3,2,1,2,2,3,1,1,3,3,2,1,2,2,3,1,1];
DeBruijn_sequence(10).vals = [1,2,2,1,3,2,3,3,1,1,2,2,1,3,2,3,3,1,1];
DeBruijn_sequence(999).vals = [1,2,1,1,3,3,2,2,3,1];
%% Find the De Bruijn sequence for each run
disp('Making contrast files...');
for r = 1:length(runNum)
    if runNum(r) == 999
        DeBruijn = DeBruijn_sequence(runNum(r)).vals;
    elseif runNum(r) > 10
        DeBruijn = DeBruijn_sequence(mod(runNum(r),10)).vals;
    else
        DeBruijn = DeBruijn_sequence(runNum(r)).vals;
    end
    names = {'1' '2' '3'};
    numTR = blockdur/TR;
    contrasts = zeros(152,3);
    for i = 1:length(DeBruijn)
        if DeBruijn(i) == 1
            contrasts(((i-1)*numTR+1:i*numTR),1) = 1;
        elseif DeBruijn(i) == 2
            contrasts(((i-1)*numTR+1:i*numTR),2) = 1;
        elseif DeBruijn(i) == 3
            contrasts(((i-1)*numTR+1:i*numTR),3) = 1;
        end
    end
    if maketext
        dlmwrite(fullfile(outputDir,['Run' num2str(runNum(r)) '_contrasts.txt']),contrasts,'delimiter','\t')
        dlmwrite(fullfile(outputDir,'baseline.txt'),contrasts(:,1),'delimiter','\t')
        dlmwrite(fullfile(outputDir,'condition_1.txt'),contrasts(:,2),'delimiter','\t')
        dlmwrite(fullfile(outputDir,'condition_2.txt'),contrasts(:,3),'delimiter','\t')
    end
end
disp('done.');