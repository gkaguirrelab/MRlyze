function [ind,Tasks,numPerTask] = break_up_matrix(matrixLength,numTasks)

% Breaks up a matrix into smaller sub-matrices, specified by 'numTasks'. 
%       Useful for running matrix calculations within a for loop.
%
%   Usage:
%   [vals,Tasks,numPerTask] = break_up_matrix(matrixLength,numTasks)
%   
%   Written by Andrew S Bock Jul 2015

%% Break up the matrix into smaller sub-matrices
[numPerTask,Tasks] = calc_tasks(matrixLength,ceil(matrixLength/numTasks));
%% Find the indices of the smaller sub-matrices
idx = [];
for i = 1:Tasks
    if isempty(idx);
        idx = [1,numPerTask(i)];
    else
        idx = [idx;[idx(end,2)+1,idx(end,2)+numPerTask(i)]];
    end
    ind{i} = idx(i,1):idx(i,2);
end