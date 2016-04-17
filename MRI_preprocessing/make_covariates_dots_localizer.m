function make_covariates_dots_localizer(session_dir,runNum,blockdur)

% Outputs covariate text files for FSL's FEAT, based on 'stimulus_display',
% and the 'dots' localizer
%
%   Usage:
%   make_covariates_dots_localizer(session_dir,runNum,blockdur)
%
%   Written by Andrew S Bock Nov 2015

%% set defaults
if ~exist('blockdur','var')
    %blockdur = 16; % duration of blocks, in seconds
    blockdur = 16.56; % SC3T, TR = 0.72
end
d = find_bold(session_dir);
out_dir = fullfile(session_dir,d{runNum});
%   The occurance of each condition is based on the De Bruijn sequence for
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
%% Get DeBruijn sequence
if runNum == 999
    DeBruijn = DeBruijn_sequence(runNum).vals;
elseif runNum > 10
    DeBruijn = DeBruijn_sequence(mod(runNum,10)).vals;
else
    DeBruijn = DeBruijn_sequence(runNum).vals;
end
%% Find timings of various frequencies
Conds = unique(DeBruijn); % Get unique frequencies
for i = 1:length(Conds)
    ct = 0;
    for j = 1:length(DeBruijn)
        if Conds(i) == DeBruijn(j)
            ct = ct + 1;
            contrasts{i}(ct,1) = blockdur*j - blockdur;
            contrasts{i}(ct,2) = blockdur;
            contrasts{i}(ct,3) = 1;
        end
    end
end
%% Write out the covariate files
for i = 1:size(contrasts,2)
    if i>1 % first condition is 0Hz flicker
        dlmwrite(fullfile(out_dir,['condition_' num2str(i-1) '.txt']),contrasts{i}(:,:),'delimiter',' ','precision','%10.5f')
    end
end