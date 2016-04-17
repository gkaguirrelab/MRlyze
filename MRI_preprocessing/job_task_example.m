function job_task_example(numTotal,numCores)

% Computes the number of tasks given the total number of interations, and
% number of cores requested.
%
%   Written by Andrew S Bock Mar 2015

%% Computer the number of tasks and number of interations per task
[numPerTask,numTasks] = calc_tasks(numTotal,numCores);
%% Create job and tasks
sched = parcluster();
job = createJob(sched);
for i = 1:numTasks
    % createTask(jobID,@some_func,numOutputs,{some_func_inputs});
    createTask(job,@foo,1,{1,numPerTask(i)});
end
submit(job);
wait(job);
out = fetchOutputs(job);
cat(2,out{:}) % concatenate all the output cells
delete(job);