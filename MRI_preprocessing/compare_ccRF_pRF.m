function compare_ccRF_pRF(xaxis,yaxis,xthresh,ythresh,area,map,areafile,eccenfile,yaxis2)

% Takes in two vectors, makes a scatter plot
%
%   Usage:
%   compare_ccRF_pRF(xaxis,yaxis,xthresh,ythresh,area,map,areafile)
%
%   Written by Andrew S Bock Nov 2014

%% set defaults
if ~exist('map','var')
    map = 'ecc';
end
if ~exist('areafile','var')
    areafile = '~/data/2014-10-29.areas-template.nii.gz';
end
if ~exist('eccenfile','var')
    eccenfile = '~/data/2014-10-29.eccen-template.nii.gz';
end
%% Load inputs
x = load_nifti(xaxis);
xt = load_nifti(xthresh);
yt = load_nifti(ythresh);
y = load_nifti(yaxis);
a = load_nifti(areafile);
e = load_nifti(eccenfile);
if exist('yaxis2','var')
    y2 = load_nifti(yaxis2);
end
%% Define ROI vertices
switch area
    case {'V1'}
        verts = find(a.vol == 1 | a.vol == -1);
    case {'V2'}
        verts = find(a.vol == 2 | a.vol == -2);
    case {'V3'}
        verts = find(a.vol == 3 | a.vol == -3);
end
%% Convert polar angle from ccRF to degrees
if strcmp(map,'polang')
    x.vol = rad2deg(x.vol)+90;
    y.vol = rad2deg(y.vol)+90;
end
%% Threshold data
% the ccRF fit is R, not R2, so square vals to match pRF
if ~isempty(strfind(xthresh,'ccRF'));
    xt.vol = xt.vol.^2;
end
if ~isempty(strfind(ythresh,'ccRF'));
    yt.vol = yt.vol.^2;
end
x.vol(xt.vol<0.15) = nan;
y.vol(yt.vol<0.15) = nan;
x.vol(e.vol>8 | e.vol<2) = nan;
y.vol(e.vol>8 | e.vol<2) = nan;
x.vol(x.vol>180) = nan;
y.vol(y.vol>180) = nan;
%% Compute correlation value
%[R,P] = corrcoef(x.vol(verts),y.vol(verts),'rows','complete'); % 'rows','complete' ignores rows with nans
%% Compute fit
if strcmp(map,'ecc')
    [slope,intercept,MSE,R2,S] = logfit(x.vol(verts),y.vol(verts),'loglog');
else
    [slope,intercept,MSE,R2,S] = logfit(x.vol(verts),y.vol(verts),'linear');
end
P = [slope,intercept];
Rinv = inv(S.R); c = (Rinv*Rinv')*S.normr^2/S.df;
se = sqrt(diag(c))';
t = P ./ se;
p = tcdf(t,S.df,'upper')*2;p=p(1);
%% Plot values using ROI vertices
figure;S = scatter(x.vol(verts),y.vol(verts),40,'filled'); hold on
%fL = [area ' ' map];
%sL = [' R2 = ' sprintf('%3.3f',R2)]; % '; P = ' sprintf('%3.3f',P(2))];
% T = title({...
%     '\fontsize{25pt}\bf{' fL '}'
%     '\fontsize{15pt}\rm{' sL '}'
% });
T = title({[area ' ' map];[' R2 = ' sprintf('%3.3f',R2) '; MSE = ' sprintf('%3.3f',MSE) '; P = ' sprintf('%3.3f',p)]});
if ~isempty(findstr('template',xaxis))
    xl = xlabel('Benson Template (degrees)');
elseif ~isempty(findstr('ccRF',xaxis))
    xl = xlabel('Subject ccRF (degrees)');
else
    xl = xlabel('Subject pRF (degrees)');
end
if ~isempty(findstr('template',yaxis))
    yl = ylabel('Benson Template (degrees)');
elseif ~isempty(findstr('ccRF',yaxis))
    yl = ylabel('Subject ccRF (degrees)');
else
    yl = ylabel('Subject pRF (degrees)');
end
set(T,'FontSize',25)
set(xl,'Interpreter','none','FontSize',20);
set(yl,'Interpreter','none','FontSize',20);
%% Compute polyfit
% X = x.vol(verts);Y = y.vol(verts);
% ind = ~isnan(X) & ~isnan(Y);
% [Po] = polyfit(X(ind),Y(ind),1);
%% Adjust plot
% Yo = polyval(Po,X(ind));plot(X(ind),Yo,'b');
if strcmp(map,'ecc')
    % Set axes
    axis('square');
    set(gca,'XScale','log');
    set(gca,'YScale','log');    
    % Set limits
    axis([0.1 8 0.1 8])
    set(gca,'xtick',0:8)
    set(gca,'ytick',0:8)
    % Manually update tick labels to have same precision.
    %   Note, they will no longer update if scale is changed
    xTick = get(gca,'xTick');
    yTick = get(gca,'yTick');
    xTickLabel = arrayfun( @(x) sprintf('%3.1f',x), xTick, 'uniformoutput', false);
    yTickLabel = arrayfun( @(x) sprintf('%3.1f',x), yTick, 'uniformoutput', false);
    set(gca,'xTickLabel',xTickLabel);
    set(gca,'yTickLabel',yTickLabel);
    % Plot fit line
    X = 0:8/100:8;
    Y = (10^intercept)*X.^(slope);plot(X,Y,'r')
    % Plot identity line
    H = refline(1,0);
    set(H,'Color',[0 0 0]);
elseif strcmp(map,'sig')
    % Set axes
    axis('square');
    axis([0 20 0 20])
    % Set limits
    
    %     set(gca,'xtick',0:180)
    %     set(gca,'ytick',0:180)
    %     % Manually update tick labels to have same precision.
    %     %   Note, they will no longer update if scale is changed
    %     xTick = get(gca,'xTick');
    %     yTick = get(gca,'yTick');
    %     xTickLabel = arrayfun( @(x) sprintf('%3.1f',x), xTick, 'uniformoutput', false);
    %     yTickLabel = arrayfun( @(x) sprintf('%3.1f',x), yTick, 'uniformoutput', false);
    %     set(gca,'xTickLabel',xTickLabel);
    %     set(gca,'yTickLabel',yTickLabel);
%     X = 0:pi/100:180;
%     Y = X*slope+intercept;plot(X,Y,'r')
    X = 0:20/100:20;
    Y = X.*(slope)+intercept;plot(X,Y,'r')
    % Plot identity line
    H = refline(1,0);
    set(H,'Color',[0 0 0]);
else
    % Set axes
    axis('square');
    % Set limits
    axis([0 180 0 180])
    %     set(gca,'xtick',0:180)
    %     set(gca,'ytick',0:180)
    %     % Manually update tick labels to have same precision.
    %     %   Note, they will no longer update if scale is changed
    %     xTick = get(gca,'xTick');
    %     yTick = get(gca,'yTick');
    %     xTickLabel = arrayfun( @(x) sprintf('%3.1f',x), xTick, 'uniformoutput', false);
    %     yTickLabel = arrayfun( @(x) sprintf('%3.1f',x), yTick, 'uniformoutput', false);
    %     set(gca,'xTickLabel',xTickLabel);
    %     set(gca,'yTickLabel',yTickLabel);
%     X = 0:pi/100:180;
%     Y = X*slope+intercept;plot(X,Y,'r')
    X = 0:pi/100:180;
    Y = x.vol(verts)*slope+intercept;plot(x.vol(verts),Y,'r'); 
    % Plot identity line
    H = refline(1,0);
    set(H,'Color',[0 0 0]);
end
if exist('yaxis2','var')
    scatter(x.vol(verts),y2.vol(verts),40,'filled','g');    
end