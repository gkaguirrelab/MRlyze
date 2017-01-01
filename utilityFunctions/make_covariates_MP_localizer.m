function make_covariates_MP_localizer(session_dir,runNum,blockdur,DeBruijn)

% Outputs covariate text files for FSL's FEAT, based on M vs P stimuli
% using the OneLight
%
%   Usage:
%   make_covariates_MP_localizer(session_dir,runNum,blockdur,DeBruijn)
%
%   Written by Andrew S Bock Nov 2015

%% set defaults
if ~exist('blockdur','var')
    blockdur = 12; % in seconds.  SC3T, TR = 0.72
end
if ~exist('DeBruijn','var')
    % 1 - Background
    % 2 - Lightflux (M stimulus - 32Hz flicker)
    % 3 - L-M (P stimulus, 2Hz flicker)
    DeBruijn = [2,1,1,3,3,2,2,3,1,2,1,1,3,3,2,2,3,1];
end
d = find_bold(session_dir);
out_dir = fullfile(session_dir,d{runNum});
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