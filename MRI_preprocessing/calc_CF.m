function [cfs] = calc_CF(srctc,trgtc,trgdists,seedSig,DoG,TR)
% Takes in a matrix of observed fmri timecourses, and outputs the
%   connective field (CF) parameters parameters that explain the most
%   variance in each observed timecourse.
%
%   Inputs:
%       srctc = TR x voxels/vertices matrix of timecourses (e.g. V2/V3
%           timecourses)
%       trgtc = TR x voxels/vertices matrix of timecourses for an ROI (e.g.
%           V1)
%       trgdists = voxel/vertices x voxel/vertices matrix of cortical
%           distances within an ROI (e.g. V1)
%       seedSig = row vector of sigma values for calculating initial
%           predicted responses (e.g. seedSig = [3 6 9]);
%
%   Usage:
%       [best] = calc_CF(srctc,trgtc,trgdists,seedSig)
%
%   Output:
%       cfs.seed = index of the best predicted timecourse for each
%           voxel/vertex
%       cfs.co = maximum correlation value for the best predicted
%           timecourse
%       cfs.sig = fit sigma value, starting with the best

%% Set the tc_shift
% postive values move predicted responses earlier in time, negative values
% move later in time
if ~exist('TR','var')
    TR = 2; % assumes 2s TR
end
tc_shift = (-1:0.5:5)/TR;
%% Convert srctc to single
srctc = single(srctc);
%%
seedX = 1:length(trgdists);
seedSig1 = seedSig{1};
seedSig2 = seedSig{2};
seedSig3 = seedSig{3};
seedSig4 = seedSig{4};
[xList] = ndgrid(seedX,seedSig1,seedSig2,seedSig3,seedSig4);
xList = xList(:);
% Create sigma list for each center location
ct = 0;
for s1 = 1:length(seedSig1);
    for s2 = 1:length(seedSig2);
        for s3 = 1:length(seedSig3);
            for s4 = 1:length(seedSig4);
                ct = ct+1;
                singleSigList(ct,:) = [seedSig1(s1),seedSig1(s1)*seedSig2(s2),seedSig3(s3),seedSig4(s4)];
            end
        end
    end
end
% Replicate this sigma list for all center locations
sigList = repmat(singleSigList,length(xList)/size(singleSigList,1),1);
% Make single (to save memory)
sigList = single(sigList);
xList = single(xList);
% Remove flat sigmas
badind = sigList(:,3) == 0 & sigList(:,4) == 0;
sigList(badind,:) = [];
xList(badind) = [];
%% Break up sigma list matrix
[sigPerTask,sigTasks] = calc_tasks(length(sigList),ceil(length(sigList)/1000));
sigidx = [];
for i = 1:sigTasks
    if isempty(sigidx);
        sigidx = [1,sigPerTask(i)];
    else
        sigidx = [sigidx;[sigidx(end,2)+1,sigidx(end,2)+sigPerTask(i)]];
    end
    sigvals{i} = sigidx(i,1):sigidx(i,2);
end
%% Initialize variables
cfs.co = zeros(size(srctc,2),1,'single');
cfs.cosig1 = zeros(size(srctc,2),1,'single');
cfs.cosig2 = zeros(size(srctc,2),1,'single');
cfs.cosig3 = zeros(size(srctc,2),1,'single');
cfs.cosig4 = zeros(size(srctc,2),1,'single');
cfs.cocenter = ones(size(srctc,2),1,'single'); % set to first index (0 will break code later)
cfs.copeakt = zeros(size(srctc,2),1,'single');
cfs.var = zeros(size(srctc,2),1,'single');
cfs.var_as_co = zeros(size(srctc,2),1,'single');
cfs.varsig1 = zeros(size(srctc,2),1,'single');
cfs.varsig2 = zeros(size(srctc,2),1,'single');
cfs.varsig3 = zeros(size(srctc,2),1,'single');
cfs.varsig4 = zeros(size(srctc,2),1,'single');
cfs.varcenter = ones(size(srctc,2),1,'single'); % set to first index (0 will break code later)
cfs.varpeakt = zeros(size(srctc,2),1,'single');
%% Create Predicted Responses
disp('Finding Connective Fields...');
for t = 1:length(tc_shift)
    disp(['Running tc_shift ' num2str(t) ' of ' num2str(length(tc_shift))]);
    tstart = clock;
    ProgressBar_parfor(tstart);
    for s = 1:sigTasks;
        [~,tmpPredResp] = calc_CF_pred(trgdists(:,xList(sigvals{s})),sigList(sigvals{s},:),trgtc,[],DoG);
        % postive values move predicted responses earlier in time, negative
        % values move later in time
        shiftPredResp = shift_tc(tmpPredResp,tc_shift(t));
        tmp_co = corr(srctc,shiftPredResp);
        tmp_co = single(tmp_co);
        % Save best co
        [cfsco,cfscoseed] = max(tmp_co,[],2);
        newcoind = cfs.co<cfsco;
        cfs.co(newcoind) = cfsco(newcoind);
        cfs.cosig1(newcoind) = sigList(sigvals{s}(cfscoseed(newcoind)),1);
        cfs.cosig2(newcoind) = sigList(sigvals{s}(cfscoseed(newcoind)),2);
        cfs.cosig3(newcoind) = sigList(sigvals{s}(cfscoseed(newcoind)),3);
        cfs.cosig4(newcoind) = sigList(sigvals{s}(cfscoseed(newcoind)),4);
        cfs.cocenter(newcoind) = xList(sigvals{s}(cfscoseed(newcoind)));
        cind = find(newcoind);
        for c = 1:length(cind);
            cfs.copeakt(cind(c)) = tc_shift(t);
        end
        % Save best var
        [cfsvar,cfsvarseed] = max(tmp_co.^2,[],2);
        newvarind = cfs.var<cfsvar;
        cfs.var(newvarind) = cfsvar(newvarind);
        cfs.varsig1(newvarind) = sigList(sigvals{s}(cfsvarseed(newvarind)),1);
        cfs.varsig2(newvarind) = sigList(sigvals{s}(cfsvarseed(newvarind)),2);
        cfs.varsig3(newvarind) = sigList(sigvals{s}(cfsvarseed(newvarind)),3);
        cfs.varsig4(newvarind) = sigList(sigvals{s}(cfsvarseed(newvarind)),4);
        cfs.varcenter(newvarind) = xList(sigvals{s}(cfsvarseed(newvarind)));
        vcind = find(newvarind);
        for vc = 1:length(vcind);
            cfs.var_as_co(vcind(vc)) = tmp_co(vcind(vc),cfsvarseed(vcind(vc)));
            cfs.varpeakt(vcind(vc)) = tc_shift(t);
        end
        ProgressBar_parfor(tstart,s,sigTasks);
    end
    ProgressBar_parfor(tstart,'clean'); % Clean up files and display total loop time
end
disp('done.');