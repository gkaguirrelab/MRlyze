function plot_template_mesh(session_dirs,template,func,cond,V2V3)
% Plots the template fits
%
%   Usage:
%   plot_template_mesh(session_dirs,template,func,cond,plot_error_bars,figsave)
%
%   Written by Andrew S Bock Oct 2015
%   Updated Jan 2016 - now only uses scatter3, shows as 2D

%% Set defaults
if ~exist('session_dirs','var')
    session_dirs = {...
        '/data/jet/abock/data/Template_Retinotopy/AEK/10012014' ...
        '/data/jet/abock/data/Template_Retinotopy/ASB/10272014' ...
        '/data/jet/abock/data/Template_Retinotopy/GKA/10152014' ...
        };
end
if ~exist('template','var')
    template = 'both';
end
if ~exist('func','var')
    func = 's5.dbrf.tf';
end
if ~exist('cond','var')
    cond = 'Movie';
end
if strcmp(template,'coarse')
    templateScale = 10;
    vals = -0.8:.1:0.8;
    tickVals = vals;
    vLabels = {'-0.8', '-0.7', '-0.6', '-0.5', '-0.4', '-0.3', '-0.2', '-0.1', '0.0', ...
        '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8'};
    GSize = 0.1; % grid sampling
elseif strcmp(template,'fine')
    templateScale = 30;
    vals = -4/30:1/30:4/30;
    tickVals = vals;
    vLabels = {'-0.133', '-0.100', '-0.066', '-0.033', '0.000', '0.033', '0.066', ...
        '0.100', '0.133'};
    GSize = 1/30; % grid sampling
elseif strcmp(template,'both')
    templateScale = 30;
    vals = -0.8:1/30:0.8;
    tickVals = -0.8:.1:0.8;
    vLabels = {'-0.8', '-0.7', '-0.6', '-0.5', '-0.4', '-0.3', '-0.2', '-0.1', '0.0', ...
        '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8'};
    GSize = 1/30; % grid sampling
end
hemis = {'lh' 'rh'};
numComps = length(session_dirs)*length(hemis);
% Plot defaults
varexpLims = [0.1 0.25];
yTicks = [0 0.05 0.10 0.15 0.2 0.25];
aLims = {[vals(1) vals(end)],[vals(1) vals(end)],varexpLims};
MSize = 10; % point marker size
LWidth = 1; % error bar line width
squareSize = 500; % areas of squares in plot
%% Load in template values
disp('Loading template values...');
if strcmp(template,'coarse') || strcmp(template,'fine')
    ct = 0;
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        for ss = 1:length(session_dirs)
            ct = ct + 1;
            session_dir = session_dirs{ss};
            if V2V3
                fitType = 'V2V3';
                tdir = fullfile(session_dir,'pRFs',template,func,cond,'V2V3');
            else
                fitType = 'V1';
                tdir = fullfile(session_dir,'pRFs',template,func,cond,'V1');
            end
            [varexp,~,~,x,y,z] = plot_template_fits(tdir,template,hemi,fitType);
            varAll(ct,:) = varexp';
            X(ct,:) = x;
            Y(ct,:) = y;
            Z(ct,:) = z;
        end
    end
elseif strcmp(template,'both')
    ct = 0;
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        for ss = 1:length(session_dirs)
            ct = ct + 1;
            session_dir = session_dirs{ss};
            % coarse
            if V2V3
                fitType = 'V2V3';
                ctdir = fullfile(session_dir,'pRFs','coarse',func,cond,'V2V3');
            else
                fitType = 'V1';
                ctdir = fullfile(session_dir,'pRFs','coarse',func,cond,'V1');
            end
            [cvarexp,~,~,cx,cy,cz] = plot_template_fits(ctdir,'coarse',hemi,fitType);
            % fine
            if V2V3
                fitType = 'V2V3';
                ftdir = fullfile(session_dir,'pRFs','fine',func,cond,'V2V3');
            else
                fitType = 'V1';
                ftdir = fullfile(session_dir,'pRFs','fine',func,cond,'V1');
            end
            [fvarexp,params,~,fx,fy,fz] = plot_template_fits(ftdir,'fine',hemi,fitType);
            fx = fx + params(1).FCx; % adjust fine center
            fy = fy + params(1).FCy; % adjust fine center
            fz = fz + params(1).psi; % adjust fine center
            % combine 'coarse' and 'fine'
            varAll(ct,:) = [cvarexp',fvarexp'];
            X(ct,:) = [cx,fx];
            Y(ct,:) = [cy,fy];
            Z(ct,:) = [cz,fz];
        end
    end
end
disp('done.');
%% Fill in subject values
sumAll = nan(numComps,length(vals),length(vals),length(vals));
sumCt = zeros(length(vals),length(vals),length(vals));
for i = 1:length(hemis)*length(session_dirs)
    for j = 1:length(X)
        tmpCoords = [X(i,j),Y(i,j),Z(i,j)];
        tmpCoords = round(tmpCoords*templateScale) + round(length(vals)/2);
        tmpVarexp = varAll(i,j);
        sumAll(i,tmpCoords(1),tmpCoords(2),tmpCoords(3)) = tmpVarexp;
        sumCt(tmpCoords(1),tmpCoords(2),tmpCoords(3)) = ...
            sumCt(tmpCoords(1),tmpCoords(2),tmpCoords(3)) + 1;
    end
end
%% Average values
avgAll = squeeze(nanmean(sumAll,1));
avgAll = avgAll(:);
%avgAll(avgAll == 0) = nan;
semAll  = squeeze(nanstd(sumAll,0,1)) ./ sqrt(sumCt);
%% Create 3D meshgrid
[sumY,sumX,sumZ] = meshgrid(vals,vals,vals);
sumX = sumX(:);
sumY = sumY(:);
sumZ = sumZ(:);
%% Plot
if V2V3
    saveName = 'V2V3';
else
    saveName = 'V1';
end
xyplane = find(sumZ == 0);
xzplane = find(sumY == 0);
yzplane = find(sumX == 0);
% XY
figure('units','normalized','position',[0 0 1 1]);
x1 = sumX(xyplane);
y1 = sumY(xyplane);
z1 = avgAll(xyplane);
e1 = semAll(xyplane);
s1 = squareSize*ones(size(x1));
ph=scatter3(x1,y1,z1,s1,z1,'s','filled');
zlim(varexpLims);
grid off
view(2)
axis square
cbh=colorbar;
colormap(hot(100));
caxis(varexpLims);
set(cbh,'YTick',yTicks);
set(gca,'XLim',[min(tickVals) max(tickVals)]);
set(gca,'XTick',tickVals);
set(gca,'XTickLabel',vLabels);
set(gca,'YLim',[min(tickVals) max(tickVals)]);
set(gca,'YTick',tickVals);
set(gca,'YTickLabel',vLabels);
xlabel('FCx','FontSize',20);
ylabel('FCy','FontSize',20);
zlabel('Variance Explained','FontSize',20);
set(gca,'FontSize',15);
savefigs('pdf',[template '_varexp_FCx_FCy_highres_' saveName]);
close all
% XZ
figure('units','normalized','position',[0 0 1 1]);
x1 = sumX(xzplane);
y1 = sumZ(xzplane);
z1 = avgAll(xzplane);
e1 = semAll(xzplane);
s1 = squareSize*ones(size(x1));
ph=scatter3(x1,y1,z1,s1,z1,'s','filled');
zlim(varexpLims);
grid off
view(2)
axis square
cbh=colorbar;
colormap(hot(100));
caxis(varexpLims);
set(cbh,'YTick',yTicks);
set(gca,'XLim',[min(tickVals) max(tickVals)]);
set(gca,'XTick',tickVals);
set(gca,'XTickLabel',vLabels);
set(gca,'YLim',[min(tickVals) max(tickVals)]);
set(gca,'YTick',tickVals);
set(gca,'YTickLabel',vLabels);
xlabel('FCx','FontSize',20);
ylabel('Psi','FontSize',20);
zlabel('Variance Explained','FontSize',20);
set(gca,'FontSize',15);
savefigs('pdf',[template '_varexp_FCx_Psi_highres_' saveName]);
close all
% YZ
figure('units','normalized','position',[0 0 1 1]);
x1 = sumY(yzplane);
y1 = sumZ(yzplane);
z1 = avgAll(yzplane);
e1 = semAll(yzplane);
s1 = squareSize*ones(size(x1));
ph=scatter3(x1,y1,z1,s1,z1,'s','filled');
zlim(varexpLims);
grid off
view(2)
axis square
cbh=colorbar;
colormap(hot(100));
caxis(varexpLims);
set(cbh,'YTick',yTicks);
set(gca,'XLim',[min(tickVals) max(tickVals)]);
set(gca,'XTick',tickVals);
set(gca,'XTickLabel',vLabels);
set(gca,'YLim',[min(tickVals) max(tickVals)]);
set(gca,'YTick',tickVals);
set(gca,'YTickLabel',vLabels);
xlabel('FCy','FontSize',20);
ylabel('Psi','FontSize',20);
zlabel('Variance Explained','FontSize',20);
set(gca,'FontSize',15);
savefigs('pdf',[template '_varexp_FCy_Psi_highres_' saveName]);
close all