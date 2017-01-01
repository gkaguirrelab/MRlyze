function calc_surface_distance(session_dir,subject_name,hemi,ROI,ROIverts,SUBJECTS_DIR)
% Calculates the distance on the cortical surface from every vertex within
%   an ROI to every other vertex within the same ROI. This is done using
%   Matlab's 'graphshortestpath' function.
%
%   Inputs:
%
%   Outputs:
%
%   Written by Andrew S Bock Apr 2015

%% Set defaults
if ~exist('hemi','var')
    hemi = 'lh';
end
if ~exist('ROI','var')
    ROI = 'prf_V1';
end
if ~exist('ROIverts','var')
    if strcmp(ROI,'V1');
        V1 = load_nifti(fullfile(session_dir,[hemi '.areas.nii.gz']));
        ROIverts = find(V1.vol<=1 & V1.vol >=-1);
    elseif strcmp(ROI,'prf_V1');
        V1 = load_nifti(fullfile(session_dir,[hemi '.areas_pRF.nii.gz']));
        ROIverts = find(V1.vol<=1 & V1.vol >=-1);
    end
end
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
anatdatadir = fullfile(SUBJECTS_DIR,subject_name);

%% Add to log
SaveLogInfo(session_dir,mfilename,session_dir,subject_name,hemi,ROI,ROIverts,SUBJECTS_DIR)

%%
%poolobj = gcp; % Gets current pool, and if no pool, creates a new one
% Load in the vertices and faces from the 'smoothwm' surface
[verts,faces] = freesurfer_read_surf(fullfile(anatdatadir,'surf',[hemi '.smoothwm']));
% find all the connected vertices
tmp = faces(:,1:2);
tmp = [tmp;faces(:,2:3)];
tmp = [tmp;faces(:,[1 3])];
% remove duplicates and inverse connections (e.g. 4-5 = 5-4)
[~,IA,~] = unique(tmp,'rows','stable');
edges = tmp(IA,:);
progBar = ProgressBar(length(edges), 'Removing inverse edges');
for i = 1:length(edges)
    tmp = edges(i,:);
    tmp = sort(tmp);
    edges(i,:) = tmp;
    if ~mod(i,1000);progBar(i);end;
end
[~,IA,~] = unique(edges,'rows','stable');
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
DG = sparse(foo1,foo2,W,length(edges),length(edges));
%DG = sparse(foo1,foo2,W,length(edges),length(edges),length(edges)*length(edges));
UG = tril(DG + DG');
% Find distances between uVertices
alldists = zeros(length(ROIverts),length(UG),'single');
%disp('Computing distances...')
%tstart = clock;
%ProgressBar_parfor(tstart,length(ROIverts));
progBar=ProgressBar(length(ROIverts),'Computing distances...');
for v = 1:length(ROIverts)
    [alldists(v,:),~,~] = graphshortestpath(UG,ROIverts(v),'directed',false,'METHOD','Dijkstra');
    %    ProgressBar_parfor(tstart,v,length(ROIverts));
    progBar(v);
end
%ProgressBar_parfor(tstart); % display total loop time
allDistances = nan(length(ROIverts),length(ROIverts));
for v = 1:length(ROIverts);
    allDistances(v,:) = alldists(v,ROIverts);
end
%delete(poolobj); % close parpool
save(fullfile(session_dir,[hemi '.' ROI '.dists.mat']),'allDistances');
disp('done.');
