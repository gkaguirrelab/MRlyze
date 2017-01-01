function make_ccRF_mat(session_dir,subject_name,template,rr,func,compute_dist,SUBJECTS_DIR,hemi,rois)

% make data structures for pRF and ccRF estimates
%
%   Usage:
%   make_ccpRF_mat(session_dir,subject,smooth,func,compute_dist,SUBJECTS_DIR,hemi)
%
%   Written by Andrew S Bock Oct 2014

%% Set default parameters
if ~exist('session_dir','var')
    error('No ''session_dir'' defined')
end
if ~exist('subject_name','var')
    error('No ''subject_name'' defined')
end
if ~exist('rr','var')
    error('No ''run'' defined')
end
if ~exist('template','var')
    template = 'anatomical'; % 'anatomical' = benson template; 'prf' = prf values
end
if ~exist('func','var')
    func = 'sdbrf.tf'; % prefix of functional files
end
if ~exist('compute_dist','var')
    compute_dist = 1;
end
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
if ~exist('hemi','var')
    hemi = {'lh' 'rh'};
end
if ~exist('rois','var')
    rois = {'V1' 'template' 'occipital' 'cortex' 'subcortical'};% defines the anatomy
end
anatdatadir = fullfile(SUBJECTS_DIR,subject_name);

%% Add to log
SaveLogInfo(session_dir,mfilename,session_dir,subject_name,template,rr,func,compute_dist,SUBJECTS_DIR,hemi,rois)

%% Find the bold run directories
d = listdir(fullfile(session_dir,'*BOLD_*'),'dirs');
if isempty(d)
    d = listdir(fullfile(session_dir,'*EPI_*'),'dirs');
end
if isempty(d)
    d = listdir(fullfile(session_dir,'RUN*'),'dirs');
end
nruns = length(d);
disp(['Session_dir = ' session_dir]);
disp(['Number of runs = ' num2str(nruns)]);
%% Create the ccRF .mat file
cd(session_dir);
Vertices = cell(length(rois),length(hemi));
uVertices = Vertices;
Vertices4plot = Vertices;
Timecourses = Vertices;
% output files
output_file = fullfile(session_dir,d{rr},[func '_' template '.mat']);
distfile = fullfile(session_dir,d{rr},[func '_' template '_allDistances.mat']); % if compute_dist, also saves a separate file with all distances
for hh = 1:length(hemi)
    % define rois
    % for each roi find vertices on the surface
    for thisroi = 1:length(rois)
        if strcmp(rois{thisroi},'V1')
            if strcmp(template,'anatomical')
                thisfile = fullfile(session_dir,[hemi{hh} '.areas.nii.gz']);
            elseif strcmp(template,'prf')
                thisfile = fullfile(session_dir,[hemi{hh} '.areas_pRF.nii.gz']);
            end
            thisbrain = load_nifti(thisfile);
            Vertices{thisroi,hh} = find(thisbrain.vol == 1 | thisbrain.vol == -1);
        elseif strcmp(rois{thisroi},'template')
            if strcmp(template,'anatomical')
                thisfile = fullfile(session_dir,[hemi{hh} '.areas.nii.gz']);
            elseif strcmp(template,'prf')
                thisfile = fullfile(session_dir,[hemi{hh} '.areas_pRF.nii.gz']);
            end
            thisbrain = load_nifti(thisfile);
            Vertices{thisroi,hh} = find(~isnan(thisbrain.vol));
        elseif strcmp(rois{thisroi},'occipital')
            thisfile = fullfile(anatdatadir,'label',[hemi{hh} '.' rois{thisroi} '.annot']);
            [~,label,colortable] = read_annotation(thisfile);
            Vertices{thisroi,hh} = find(label==colortable.table(10));
        elseif strcmp(rois{thisroi},'subcortical')
            thisfile = fullfile(session_dir,d{rr},[func '.nii.gz']);
            tmp = load_nifti(thisfile);
            dims = size(tmp.vol);
            % run pRF analysis for slices in the middle half of the volume
            xslices = round(dims(1)/4):round(3*dims(1)/4);
            yslices = round(3*dims(2)/6):round(4*dims(2)/6);
            zslices = round(dims(3)/6):round(3*dims(3)/6);
            ind = zeros([dims(1),dims(2),dims(3)]);
            ind(xslices,yslices,zslices) = 1;
            Vertices{thisroi,hh} = find(ind);
        elseif strcmp(rois{thisroi},'cortex')
            % just need to load in a file
            thisfile = fullfile(session_dir,[hemi{hh} '.areas.nii.gz']);
            thisbrain = load_nifti(thisfile);
            % indexes all vertices on the surface
            Vertices{thisroi,hh} = (1:length(thisbrain.vol))';
        end
    end
    %% load timecourses for all vertices
    times_roi = cell(length(rois),2);
    for thisroi = 1:length(rois)
        if strcmp(rois{thisroi},'subcortical')
            thisfile = fullfile(session_dir,d{rr},[func '.nii.gz']);
            disp([rois{thisroi} ' timecourses from ' thisfile])
            tmp = load_nifti(thisfile);
            allbrain = reshape(tmp.vol,size(tmp.vol,1)*size(tmp.vol,2)*size(tmp.vol,3),size(tmp.vol,4));
            times_roi{thisroi,hh} = allbrain(Vertices{thisroi,hh},:)';
        else
            thisfile = fullfile(session_dir,d{rr},[func '_surf.' hemi{hh} '.nii.gz']);
            disp(['Loading ' rois{thisroi} ' timecourses from ' thisfile])
            tmp = load_nifti(thisfile);
            allsurface = tmp.vol;
            times_roi{thisroi,hh} = squeeze(allsurface(Vertices{thisroi,hh},:,:,:))';
        end
    end
    %% find and save unique timecourses
    for thisroi = 1:length(rois)
        % mean and std of all timecourses (shortcut to use unique below)
        allvx = mean(times_roi{thisroi,hh}) + std(times_roi{thisroi,hh});
        % find the unique timecourses
        [uniquetc,all2unique] = unique(allvx);
        % save the unique timecourses
        Timecourses{thisroi,hh} = times_roi{thisroi,hh}(:,all2unique,:);
        % save the corresponding vertices
        uVertices{thisroi,hh} = Vertices{thisroi,hh}(all2unique);
        % for each 'unique vertex' save all vertices that have identical timecourses (for plotting)
        twins = [];
        progBar = ProgressBar(length(uniquetc),['Finding unique timecourses for ' rois{thisroi}]);
        for i = 1:length(uniquetc)
            % identify twins
            twinvertices = find(allvx == uniquetc(i));
            % also save for later
            twins = [twins; twinvertices'];
            % set up struct
            Vertices4plot{thisroi,hh}.voxel(i).uvox = uVertices{thisroi,hh}(i); %
            Vertices4plot{thisroi,hh}.voxel(i).vertices = Vertices{thisroi,hh}(twinvertices); %
            Vertices4plot{thisroi,hh}.voxel(i).n = length(twinvertices);
            if ~mod(i,10000);progBar(i);end
        end
        % note: the first one contains all flat timecourses
    end
    %% compute distance
    if compute_dist
        poolobj = gcp; % Gets current pool, and if no pool, creates a new one
        clear tmp edges verts faces DG UG
        allDistances{1,hh}=nan(length(uVertices{1,hh}),length(uVertices{1,hh}));
        % Load in the vertices and faces from the 'smoothwm' surface
        [verts,faces] = freesurfer_read_surf(fullfile(anatdatadir,'surf',[hemi{hh} '.smoothwm']));
        % find all the connected vertices
        tmp = faces(:,1:2);
        tmp = [tmp;faces(:,2:3)];
        tmp = [tmp;faces(:,[1 3])];
        % remove duplicates and inverse connections (e.g. 4-5 = 5-4)
        [C,IA,IC] = unique(tmp,'rows','stable');
        edges = tmp(IA,:);
        progBar = ProgressBar(length(edges), 'Removing inverse edges');
        for i = 1:length(edges)
            tmp = edges(i,:);
            tmp = sort(tmp);
            edges(i,:) = tmp;
            if ~mod(i,1000);progBar(i);end;
        end
        [C,IA,IC] = unique(edges,'rows','stable');
        edges = edges(IA,:);
        % Find the distance from every vertex to it's neighbors
        W = zeros(1,length(edges));
        progBar = ProgressBar(length(edges), 'Finding the distance from every vertex to it''s neighbors');
        for e = 1:length(edges)
            p = verts(edges(e,1),:);
            q = verts(edges(e,2),:);
            W(e) = sqrt((p(1)-q(1))^2+(p(2)-q(2))^2+(p(3)-q(3))^2);
            if ~mod(e,1000);progBar(e);end
        end
        % Create sparse matrix for distance analysis
        foo1 = edges(:,1)';foo2 = edges(:,2)';
        DG = sparse(foo1,foo2,W,length(edges),length(edges),length(edges)*length(edges));
        UG = tril(DG + DG');
        % Find distances between uVertices
        alldists = zeros(length(uVertices{1,hh}),length(UG),'single');
        disp('Computing distances...')
        tstart = clock;
        ProgressBar_parfor(tstart);
        %progBar = ProgressBar(length(uVertices{1,hh}),'computing distances...');
        parfor v = 1:length(uVertices{1,hh})
            [alldists(v,:),~,~] = graphshortestpath(UG,uVertices{1,hh}(v),'directed',false,'METHOD','Dijkstra');
            ProgressBar_parfor(tstart,v,length(uVertices{1,hh}));
            %if ~mod(v,5);progBar(v);end
        end
        ProgressBar_parfor(tstart,'clean'); % Clean up files and display total loop time
        delete(poolobj); % close parpool
        if compute_dist
            for v = 1:length(uVertices{1,hh})
                allDistances{1,hh}(v,:) = alldists(v,uVertices{1,hh});
            end
        end
    end
    %     % load previously saved file
    if ~compute_dist
        load(distfile,'allDistances')
    end
    %% load templates
    if strcmp(template,'anatomical')
        polfile = fullfile(session_dir,[hemi{hh} '.pol.nii.gz']);
    elseif strcmp(template,'prf')
        polfile = fullfile(session_dir,[hemi{hh} '.pol_pRF.nii.gz']);
    end
    tmp = load_nifti(polfile);
    allsurface = tmp.vol;
    Pol{1,hh} = allsurface(uVertices{1,hh});
    if strcmp(template,'anatomical')
        eccfile = fullfile(session_dir,[hemi{hh} '.ecc.nii.gz']);
    elseif strcmp(template,'prf')
        eccfile = fullfile(session_dir,[hemi{hh} '.ecc_pRF.nii.gz']);
    end
    tmp = load_nifti(eccfile);
    allsurface = tmp.vol;
    Ecc{1,hh} = allsurface(uVertices{1,hh});
end
%% Show variables
uVertices
Vertices
Vertices4plot
Timecourses
Pol
Ecc
allDistances
%%
disp('saving data file...');
save(output_file,'uVertices','Vertices','Vertices4plot','Timecourses','Pol','Ecc','allDistances','-v7.3')
if compute_dist
    save(distfile,'allDistances','-v7.3')
end
