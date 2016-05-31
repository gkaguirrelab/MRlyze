function [hcamlight] = surface_plot(map_type,data,subject_name,hemi,surface,thresh,trans,polar_map,show_curv,view_angle,light_angle)
%   Visualize a surface overlay in matlab
%
%   Usage:
%   surface_plot(map_type,data,subject_name,hemi,surface,thresh,trans,polar_map,view_angle,light_angle)
%
%   map_types:
%       'ecc' 'pol' 'sig' 'co' 'var' 'areas' 'zstat'
%   data = full path to data to overlay on surface or data array
%       (note, must be in the same space as 'surface')
%   subject_name = subject surface to plot (e.g. 'fsaverage')
%   hemi = hemisphere ('lh', 'rh');
%   surface = 'inflated' % 'pial' or any surface in freesurfer
%   trans = transparancy value (0.65 - default)
%   thresh = logical vector of the same length as data file, 0s for
%   vertices not to show (e.g. less than some threshold).
%   polar_map = portion of visual field. 'hemi' or 'full' ('hemi' - default);
%   view_angle = view angle (default: [45,-10] for 'lh' and [-45,-10] for 'rh');
%   light_angle = light angle (default: [45,-10] for 'lh' and [-45,-10] for 'rh');
%
%   example: surface_plot('ecc','/home/andrew/SUBJECTS/RESTING/data/stim_lh_occipital_resting_ecc_Templ_rf.tf_avgsurf.nii.gz','fsaverage','lh',0.75)
%
%   written by Andrew Bock 2013

%% Set up defaults
SUBJECTS_DIR=getenv('SUBJECTS_DIR');
if ~exist('subject_name','var')
    subject_name = 'fsaverage_sym';
end
if ~exist('hemi','var')
    hemi = 'lh';
end
if ~exist('surface','var')
    surface = 'inflated';
end
if ~exist('trans','var')
    trans = 0.85;
end
if ~exist('polar_map','var')
    polar_map = 'full'; %'hemi'
end
if ~exist('view_angle','var')
    if strcmp(hemi,'lh')
        view_angle = [30,-10];
    elseif strcmp(hemi,'rh')
        view_angle = [-30,-10];
    end
end
if ~exist('show_curv','var')
    show_curv = 1;
end
if ~exist('light_angle','var')
    if strcmp(hemi,'lh')
        light_angle = [30,-10];
    elseif strcmp(hemi,'rh')
        light_angle = [-30,-10];
    end
end
%% Define color map
switch map_type
    case 'ecc'
        mapres=[0 15 200];
        %mycolormap = blue_red_yellow;
        mycolormap = make_ecc_colormap(mapres);
    case 'pol'
        mapres=[-pi pi 200];
        if strcmp(polar_map,'hemi')
            mycolormap = make_polar_colormap_hemi(mapres);
        else
            mycolormap = make_polar_colormap(mapres);
        end
    case 'sig'
        mycolormap = flipud(jet(200));
        mapres=[0 10 200];
    case 'peak'
        mycolormap = flipud(jet(200));
        mapres=[3 10 200];
    case 'co'
        %         mycolormap = jet(200);
        mapres=[-1 1 200];
        mycolormap = make_rbco_colormap(mapres);
    case 'rbco'
        mapres=[-1 1 200];
        mycolormap = make_rbco_colormap(mapres);
    case 'var'
        mycolormap = hot(200);
        mapres=[0 1 200];
    case 'invvar'
        mycolormap = flipud(hot(200));
        mapres=[0 1 200];
    case 'DoG'
        mycolormap = jet(200);
        mapres=[-.5 1 200];
    case 'areas'
        mycolormap = [
            .75 .75 .75
            1 1 1
            .75 .75 .75
            1 1 1
            .75 .75 .75
            1 1 1
            .75 .75 .75
            1 1 1
            .75 .75 .75
            1 1 1
            ];
        mapres=[-5 5 length(mycolormap)];
    case 'blueareas'
        mycolormap = [
            0 0 .5
            0 .33 1
            0 .75 1
            0 .75 1
            0 .33 1
            0 0 .5
            ];
        mapres=[-3 3 length(mycolormap)];
    case 'ROI'
        mycolormap = [
            1 1 0
            ];
        mapres=[1 1 length(mycolormap)];
    case 'zstat'
        %         mycolormap = jet(200);
        mapres=[-10 10 200];
        mycolormap = make_rbco_colormap(mapres);
    case 'Fstat'
        mycolormap = hot(200);
        mapres=[0 100 200];
    case 'eccdiff'
        mycolormap = jet(200);
        mapres=[0 10 200];
    case 'poldiff'
        mycolormap = jet(200);
        mapres=[0 2*pi 200];
    case 'residual'
        mycolormap = hot(200);
        mapres=[0 200 200];
    case 'diffecc'
        mapres=[0 10 200];
        %mycolormap = blue_red_yellow;
        mycolormap = make_ecc_colormap(mapres);
    case 'diffpol'
        mapres=[0 pi 200];
        %mycolormap = blue_red_yellow;
        mycolormap = make_ecc_colormap(mapres);
    case 'logp'
        mycolormap = flipud(hot(200));
        mapres=[-5 0 200];
    case 'other'
        mycolormap = make_colormap([1 1 0;1 0 0;0 1 0;0 0 1;1 0 1]);
        mapres=[0 21 size(mycolormap,1)];
    otherwise
        disp('map not recognized');
end
if strcmp(map_type,'image')
    % elseif strcmp(map_type,'areas') || strcmp(map_type,'ROI')
    %     myvec = linspace(mapres(1),mapres(2),length(mycolormap));
else
    myvec = linspace(mapres(1),mapres(2),size(mycolormap,1));
end
%mycolormap(1,:) = [1 1 1]*.8; % set the first value of the colormap to gray (either exactly 0, which would be bad luck, or NaN)
% set the first value of the colormap to gray (either exactly 0, which would be bad luck, or NaN)
%% Load surface and curvature files
[vert,face] = freesurfer_read_surf(fullfile(SUBJECTS_DIR,subject_name,'surf',[hemi '.' surface]));
if strcmp(surface(1:3),'0.1') || strcmp(surface(1:3),'0.0')
    [curv,~] = freesurfer_read_curv(fullfile(SUBJECTS_DIR,subject_name,'surf',[hemi '.' surface '.curv']));
else
    [curv,~] = freesurfer_read_curv(fullfile(SUBJECTS_DIR,subject_name,'surf',[hemi '.curv']));
end
mycurv = -curv;
if show_curv
    ind.sulci = mycurv<0;
    ind.gyri = mycurv>0;
    ind.medial = mycurv==0;
    %     mycurv(ind.sulci) = .5;
    %     mycurv(ind.gyri) = 0.85;
    mycurv(ind.sulci) = .8;
    mycurv(ind.gyri) = 0.9;
    mycurv(ind.medial) = 0.7;
else
    mycurv = 0.85 * ones(size(mycurv));
end
%mycurv(:) = 0.375; % for making entire cortex gray
cmap_curv = repmat(mycurv,1,3);
% put all into a patch structure
brain.vertices = vert;
brain.faces = face;
brain.facevertexcdata = cmap_curv;
%% Load map data file
if ischar(data)
    if exist(data,'file')
        [~,~,ext] = fileparts(data);
        if strcmp(ext,'.gz')
            tmp = load_nifti(data);
            srf = tmp.vol;
        elseif strcmp(ext,'.mgh') || strcmp(ext,'.mgz')
            srf = load_mgh(data);
        end
    else
        error('data input does not correspond to an existing file');
    end
else
    if isequal(size(vert,1),size(data,1))
        srf = double(data);
        data = 'input array'; % Specify a string to display as the title of the figure
    else
        error(['The size of ''data'' does not correspond to ' subject_name ' surface size']);
    end
end
%% Threshold data
%srf(srf<thresh) = nan;
% need to fix this so it thresholds based on another file
%% Color vertices
cmap_vals = zeros(size(cmap_curv))+0.5;
alpha_vals = zeros(size(cmap_curv,1),1);
if strcmp(map_type,'image')
    cmap_vals = squeeze(srf);
    for i = 1:length(srf)
        if isnan(cmap_vals(i,1));
            cmap_vals(i,:) = [.8 .8 .8];
        end
    end
else
    for i = 1:length(srf)
        % Find the closest color value to the srf(i) value
        [~,ind] = min(abs(myvec-srf(i)));
        if isnan(srf(i))
            col4thisvox = [.8 .8 .8]; % set nan to gray
        else
            col4thisvox = mycolormap(ind,:);
        end
        cmap_vals(i,:) = col4thisvox;
    end
end
%% Set transparency
if strcmp(map_type,'image')
    alpha_vals(~isnan(srf(:,1))) = trans;
else
    alpha_vals(~isnan(srf)) = trans;
end
if exist('thresh','var')
    alpha_vals(~thresh) = 0;
end

% if strcmp(map_type,'pol')
%     if strcmp(polar_map,'hemi')
%         alpha_vals(srf < -pi/2 | srf > pi/2) = 0;
%     end
% end
%% Make figure
figure('units','normalized','position',[0 0 1 1]);
if strcmp(map_type,'image')
else
    % make legend first, so that brain surface is active last
    subplot(4,4,4);
    if strcmp(map_type,'ecc') || strcmp(map_type,'pol')
        x = linspace(-1,+1,mapres(3)); [xx,yy] = meshgrid(x);
        [th,r] = cart2pol(xx,yy);
        if strcmp(map_type,'ecc')
            img = r;
            img(r>1) = nan;
        elseif strcmp(map_type,'pol')
            img = th;
            if strcmp(polar_map,'hemi')
                img(r>1) = -3;
            else
                img(r>1) = 3.2;
            end
            %img(r>1) = nan;
        else
            img = yy - min(yy(:));
        end
        img = img/nanmax(img(:))*mapres(2);
        img(isnan(img)) = mapres(2);
        imagesc(x,x,img)
    end
    if ~strcmp(map_type,'pol') && ~strcmp(map_type,'ROI') && ~strcmp(map_type,'image')
        cb = colorbar;
        caxis([mapres(1) mapres(2)]);
    end
    axis tight equal off
    % if strcmp(map_type,'pol')
    %     legendmap = [.8 .8 .8;.8 .8 .8;mycolormap];
    % else
    if strcmp(map_type,'areas') || strcmp(map_type,'blueareas') ...
            || strcmp(map_type,'co') || strcmp(map_type,'rbco')
        legendmap = mycolormap;
    else
        if strcmp(polar_map,'hemi') && strcmp(map_type,'pol')
            legendmap = mycolormap;
        else
            legendmap = [mycolormap;.8 .8 .8];
        end
    end
    % end
    if strcmp(map_type,'image')
    else
        colormap(legendmap)
    end
    % Flip map if polar angle (e.g. upper visal field - ventral surface)
    % if strcmp(map_type,'pol')
    %     set(cb,'YDir','reverse')
    % else
    %     set(cb,'YDir','normal')
    % end
end
%% Plot brain and surface map
smp = brain;
smp.facevertexcdata = cmap_vals;
set(gcf,'name',data);
subplot(1,3,2); hold on
hbrain = patch(brain,'EdgeColor','none','facecolor','interp','FaceAlpha',1);
hmap = patch(smp,'EdgeColor','none','facecolor','interp','FaceAlpha','flat'...
    ,'FaceVertexAlphaData',alpha_vals,'AlphaDataMapping','none');
daspect([1 1 1]);
% Camera settings
cameratoolbar;
camproj perspective; % orthographic; perspective
lighting phong; % flat; gouraud; phong
material dull; % shiny; metal; dull
view(view_angle(1),view_angle(2)); % ASB/AEK (45,-10); GKA (55,-5)
%lightangle(light_angle(1),light_angle(2));
hcamlight = camlight('headlight');
axis tight off
zoom(2)
% note - to delete light, type 'delete(findall(gcf,'Type','light'))'