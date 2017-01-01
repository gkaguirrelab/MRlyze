function [F,df] = F_test(mat)

% Computes an F-statistic, given an input matrix.  Rows correspond to
% subjects/runs, columns to conditions.
%
%   Usage:
%   [F,df] = F_test(mat)
%
%   Written by Andrew S Bock Aug 2015

%% Pull out values
means.cond = squeeze(mean(mat)); % condition means
means.sub = squeeze(mean(mat,2)); % subject means
tmpmat = reshape(mat,size(mat,1)*size(mat,2),size(mat,3));
means.tot = mean(tmpmat,1); % grand mean
numConds = size(mat,2);
numSubs = size(mat,1);
% sum of squares conditions
for i = 1:numConds
    tmpc(i,:) = (means.cond(i,:) - means.tot).^2;
end
SSa = numSubs * sum(tmpc);
% sum of squares subjects
for i = 1:numSubs
    tmps(i,:) = (means.sub(i,:) - means.tot).^2;
end
SSs = numConds * sum(tmps);  
% sum of squares total
for i = 1:numConds*numSubs
    tmpt(i,:) = (tmpmat(i,:) - means.tot).^2;
end
SSt = sum(tmpt);
% Sum of squares error
SSas = SSt - (SSa + SSs);
%% Calculate F-stat
MSa = SSa / (numConds-1);
MSas = SSas / ((numConds-1) * (numSubs-1));
F = MSa ./ MSas;
df = [(numConds-1) (numSubs-1)];
