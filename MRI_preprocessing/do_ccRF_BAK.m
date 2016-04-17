function do_ccRF(session_dir,template,rr,func,roi,hemi,dataDir,toolboxDir)

%   Usage:
%   do_ccRF(session_dir,template,rr,func,roi,hemi,dataDir,toolboxDir)
%
%   session_dir  - contains all bold runs
%   rr - individual run to analyze (e.g. 1)
%   func - sdbrf.tf <default>
%   roi  - 3 <default>
%
% create 'true ccRFs' (generate responses, set options, etc) to recover with pRF_Estimate
% 'connection fields analysis' i.e. input = V1 activity,
% two methods: V1 activity as function of template ecc and polang;
% or as function of cortical distance (Haak 2013)
%
% PB 03/2013
% modified AB 04/2013

%% Define initial variables
if ~exist('session_dir','var')
    error('No ''session_dir'' defined')
end
if ~exist('template','var')
    error('No ''template'' defined'); % 'anatomical' or 'prf'
end
if ~exist('rr','var')
    error('No ''run'' defined')
end
if ~exist('func','var')
    func = 'sdbrf.tf'; % prefix of functional files
end
if ~exist('roi','var')
    roi = 3;% 1 - V1; 2 - V1-V3 template; 3 - occipital; 4 - cortex; 5 - subcortical
end
if ~exist('hemi','var');
    hemi = {'lh','rh'};
end
if ~exist('dataDir','var')
    dataDir = '/Users/abock/data/Retinotopy/';
end
if ~exist('toolboxPath','var')
    toolboxDir = '/Users/Shared';
end
%% Add to log
SaveLogInfo(session_dir,mfilename,session_dir,template,rr,func,roi,hemi,dataDir,toolboxDir)

%% Add files/folders to path
toolboxPath = genpath(toolboxDir);
addpath(toolboxPath);
dataPath = genpath(dataDir);
addpath(dataPath);
%% Find the bold runs in the session_dir
d = listdir(fullfile(session_dir,'*bold_*'),'dirs');
% defines the anatomy
rois = {'V1' 'V1-V3_template' 'occipital' 'cortex' 'subcortical'};% 'brain'};%{'V1' 'V2' 'V3' 'brain'};
disp(rois{roi});
%analysis_type = {'Templ' 'Dist'};
analysis_type = {'ccRF'}; %'Templ' or 'Dist'
sig_type = 0; % 0 for linear, 1 for logrithmic
%for rr = 1:nruns
if ischar(rr)
    rr = str2double(rr);
end
datafile = fullfile(session_dir,d{rr},[func '_' template '.mat']);
disp(datafile);
%input('Is the subject, datafile, and corrthr sigma type correct?')
for h =1:length(hemi)%
    for aa=1
        clear pRFs HDR opt misc seeds tRF vtc sim s Ecc Pol Vertices4plot uVertices
        prfs_dir = fullfile(session_dir,'prfs',[func '_' template]);
        system(['mkdir ' fullfile(session_dir,'prfs')]);
        system(['mkdir ' prfs_dir]);
        % set path to output file (where prfs are saved)
        output_filename = fullfile(prfs_dir,[hemi{h} '_' func ...
            '_' template '_' rois{roi} '_run' num2str(rr)  '.mat']);
        % load timecourses (etc) data -- created by make_data
        disp('Loading data file');
        load(datafile)
        
        sim = 0; % 1=simulated or 0=real data
        
        % selected runs
        runs = '1:size(Timecourses{1,1},3)'; % '1'; %'1:2'; %
        if sim==1, runs = '1'; end
        
        % set optional inputs to pRF_Estimate [type help pRF_Estimate]
        opt.parallel = 1; % 0 if running on the cluster, 1 if local
        opt.corrthr = 0.1;
        opt.freelist = {'sig'};
        opt.normalize = 1; % Convert to percent signal change
        %opt.corrthr = Inf; % if Inf only finds best seed, no nonlinear fit
        % opt.gof = 'mse';
        %opt.freelist = {'sig'}; % to fit only sigma values
        
        
        %% Convert to percent signal change
        if opt.normalize
            fprintf('Converting timecourses to percent signal change...\n')
            zeroct = 0;
            for r = [1,roi]
                tmp = Timecourses{r,h};
                for voxNum = 1:size(tmp,2)
                    ObsResp = tmp(:,voxNum);
                    if mean(ObsResp) < 1;
                        %                             zeroct = zeroct + 1;
                        %                             if zeroct >10; % more than 10 voxels with zero mean
                        %                                 disp(['hemi = ' num2str(hh)])
                        %                                 disp(['ROI = ' num2str(r)]);
                        %                                 disp(['mean = ' num2str(mean(ObsResp))]);
                        %                                 disp(['voxNum = ' num2str(voxNum)]);
                        %                                 error('Voxel mean is near zero, possibly converted to percent signal change already')
                        %                             end
                    else
                        dc = ones(size(ObsResp,1),1)*mean(ObsResp);
                        newObsResp = ((ObsResp./dc) - 1) .*100;
                        tmp(:,voxNum) = newObsResp;
                    end
                end
                Timecourses{r,h} = tmp;
            end
        end
        %% 'stimulus' (i.e. V1) space definition
        
        % s -> V1 data
        %     .x,.y     x,y from eccentricity & polar angle
        %     .nx       length x
        %     .ny   	1
        %     .dx,.dy   1,1
        %     .id       surface vertex identifier (for plotting) -- this is an extra-field, not present for in stimbased analysis
        %     .frames	timecourses [1:nVx_V1, 1, nTRs, runs]
        
        s.id = uVertices{1,h};
        th = Pol{1,h};
        r = Ecc{1,h};
        s.ecc = r;
        s.pol = th;
        [s.x,s.y] = pol2cart(s.pol,s.ecc);
        s.nx = length( s.x ) ;
        s.ny = 1 ;
        s.dx = 1;
        s.dy = 1;
        
        %% 'stimulus' (i.e. V1) timecourse
        theseruns = eval(runs);
        mymat = Timecourses{1,h}; %[time,vx,runs]
        for r = 1:size(Timecourses{1,h},3)
            for t = 1:size(mymat,1)
                s.frames(:,1,t,r) = mymat(t,:,r); %[vx,1,time,runs]
            end
        end
        theseruns = eval(runs);
        s.frames = s.frames(:,:,:,theseruns);
        stim.numTRs = size(s.frames,3);
        stim.secPerTR = 1; % never used unless there is convolution
        HDR.dt = 1;
        HDR.function = 1; % this way there is no convolution going on
        opt.hdr = HDR;
        
        % unwrap stimulus timecourse [time,space1d,nruns]
        UnwrappedConvStim = convUnwrap(s,HDR,stim);
        
        %% seeds (and their predicted timecourses)
        disp('Creating predicted timecourses');
        switch analysis_type{aa}
            case 'ccRF'
                seedX = 1:s.nx;
                    seedSig = [3 6 9];
                    %seedSig = 2;
                %seedSig = 0.5;
                [xList,sigList] = ndgrid(seedX,seedSig);
                xList = xList(:);
                sigList = sigList(:);
                nseeds = length(xList);
                RF.center = [0 0];
                s_tmp.y = 0;
                seeds.predmat = NaN(size(UnwrappedConvStim,1),size(UnwrappedConvStim,3),nseeds);
                for i = 1:nseeds
                    s_tmp.x = allDistances{1,h}(xList(i),:)';
                    RF.sig = sigList(i);
                    [~,seeds.predmat(:,:,i)] = fitPRF_AB(RF,UnwrappedConvStim,[],s_tmp,opt);
                    seeds.rfs(i).center = [xList(i), 0];
                    seeds.rfs(i).sig = [sigList(i) NaN];
                end
            case 'Templ'
                p = [0.0401 0.1331] ; % cortical magnification parameters,
                % roughly according to cortical magnification, based  on
                %  Dougherty and Wandell 2003 JOV
                if sig_type
                    seedSig = exp(linspace(log(5), log(15), 5)); % seed sigma scales with eccentricity
                else
                    seedSig = linspace(1,5,5); % sigma linear scale
                end
                step = 1; % one every x values
                seedX = s.x(1:step:end);
                [xList,sigList] = meshgrid(seedX,seedSig);
                seedY = s.y(1:step:end);
                [yList,~] = ndgrid(seedY,seedSig);
                xList = xList(:);
                yList = yList(:);
                sigList = sigList(:);
                for i = 1:length(xList)
                    seeds.rfs(i).center = [xList(i),yList(i)];
                    seeds.rfs(i).ecc = norm(seeds.rfs(i).center);
                    if sig_type
                        seeds.rfs(i).sig = [sigList(i) * p(1) * seeds.rfs(i).ecc + p(2) NaN]; % sigma scales with eccentricity
                    else
                        seeds.rfs(i).sig = [sigList(i) NaN]; % sigma linear scale
                    end
                    % line from Ione's code
                    %     Loc(hh).sigS=double(Loc(hh).sig.*((p(1).*(Loc(hh).ecc.*180/pi))+p(2)));
                end
        end
        %% ROI activity
        myROIvx = uVertices{roi,h};
        ROIindices = 1:length(myROIvx);
        nVx = size(Timecourses{roi,h},2);
        % load timecourses
        for vx = 1:nVx
            % save vx index
            vtc(vx).id = myROIvx(vx);
            % save the timecourse that is passed to pRF_Estimate
            vtc(vx).tc = squeeze(Timecourses{roi,h}(:,ROIindices(vx),:));
            % select only the runs to be analyzed
            theseruns = eval(runs);
            vtc(vx).tc = vtc(vx).tc(:,theseruns);
        end
        tRF = [];
        %% call pRF_estimate (and open parallel pool, if an option)     
        % set up opt
        opt.allDistances = allDistances;
        opt.xList = xList;
        opt.hh = h;
        
        disp('calling pRF_estimate function...')
        if ~exist('seeds','var')
            seeds = [];
        end
        
        %             if opt.parallel
        %                 %ClusterInfo.setQueueName('Default');
        %                 job = createMatlabPoolJob();
        %                 task = createTask(job, @pRF_Estimate_AB, 6, {'vtc','stim','s','opt','seeds'});
        %                 job.MinimumNumberOfWorkers = 12;
        %                 job.MaximumNumberOfWorkers = 12;
        %                 tic; submit(job); wait(job); toc
        %                 results = getAllOutputArguments(job);
        %                 save(output_filename, 'results','tRF', 'vtc','sim','s','Ecc','Pol','Vertices4plot','uVertices','-v7.3');
        %             end
        
        %j = batch('pRF_Estimate_AB','matlabpool',10,'CaptureDiary',true);
        
        %         for voxNum = 1:length(vtc)
        %             % Set up lock, done, and error files
        %             opt.lockfile(voxNum) = fullfile(session_dir,'mgo_files',['lock_' num2str(voxNum)]);
        %             opt.donefile(voxNum) = fullfile(session_dir,'mgo_files',['done_' num2str(voxNum)]);
        %             opt.errorfile(voxNum) = fullfile(session_dir,'mgo_files',['error_' num2str(voxNum)]);
        %         end
        
        if opt.parallel
            %if strcmp(analysis_type,'ccRF')
            [pRFs,HDR,opt,misc,seeds,corrvals] = pRF_Estimate_AB(vtc,stim,s,opt,seeds);
            %elseif strcmp(analysis_type,'Templ')
            %[pRFs HDR opt misc seeds] = pRF_Estimate(vtc,stim,s,opt,seeds);
            %end
            cd(session_dir)
            save(output_filename, 'pRFs', 'HDR', 'opt', 'misc', 'seeds', 'corrvals', ...
                'tRF', 'vtc','sim','s','Ecc','Pol','Vertices4plot','uVertices','-v7.3');
        end
        
    end
end