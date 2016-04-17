function [integerPerTask, numTasks] = calc_tasks(intVal, numTasks)
%PCTDEMO_HELPER_SPLIT_SCALAR Divides a non-negative integer into a sum of
%smaller non-negative integers.
%   [integerPerTask, numTasks] = PCTDEMO_HELPER_SPLIT_SCALAR(intVal, numTasks)
%   assigns a vector of length min(numTasks, intVal) to integerPerTask.
%   The sum of that vector is intVal.
%   The value of numTasks returned equals min(numTasks, numIntVal).
%
%   The input arguments must be integers greater than or equal to zero.  If
%   intVal is greater than zero, numTasks must be greater than zero.
%
%   The function is useful when dividing a Monte-Carlo simulation that is
%   repeated intVal times into numTasks tasks.  In that case, task i should
%   perform integerPerTask(i) simulations.
%
%   See also PCTDEMO_HELPER_SPLIT_VECTOR

%   Copyright 2007-2012 The MathWorks, Inc.

%% Validate the input arguments.
narginchk(2, 2);
tc = TypeChecker();
if ~(tc.isIntegerScalar(intVal, 0, Inf) ...
        && tc.isIntegerScalar(numTasks, 0, Inf))
    error('pctexample:splitscalar:SplitScalarInputsMustBePositive', ...
        'Input arguments must be non-negative integers');
end
if (intVal > 0 && numTasks == 0)
    error('pctexample:splitscalar:SplitScalarInvalidNumTasks', ...
        ['Number of tasks must be greater than 0 if the scalar is '...
        'greater than 0']);
end
% Input arguments have been validated.
if (intVal < numTasks)
    numTasks = intVal;
end
if (intVal == 0)
    integerPerTask = [];
    return;
end
% At this point, both intVal and numTasks are strictly positive integers.
split = fix(intVal / numTasks);
remainder = intVal - numTasks * split;
integerPerTask = zeros(numTasks, 1);
integerPerTask(:) = split;
integerPerTask(1:remainder) = integerPerTask(1:remainder) + 1;
end % End of pctdemo_helper_split_scalar
