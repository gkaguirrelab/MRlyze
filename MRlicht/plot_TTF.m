function plot_TTF(means,SEMs,ROI,condName,hemi,func,xlims,ylims,xTick,xLabels)

% Plots TTFs, which result from 'make_TTF'
%
%   Usage:
%   plot_TTF(means,SEMs,ROI,condName,hemi,func,xlims,ylims,xTick,xLabels)
%
%   Written by Andrew S Bock Sep 2015
%   March 2016: updated to include significant info in title. GF 
%% set defaults
if ~exist('condName','var')
    condName = '';
end
if ~exist('hemi','var')
    hemi = '';
end
if ~exist('func','var')
    func = '';
end
if ~exist('ROI','var')
    ROI = 'LGN';
end
if ~exist('xlims','var')
    xlims = [0 7];
end
if ~exist('ylims','var')
    ylims = [-1 2.5];
end
if ~exist('xTick','var')
    xTick = [0 1 2 3 4 5 6 7];
end
if ~exist('xLabels','var')
    xLabels = {'0Hz','2Hz','4Hz','8Hz','16Hz','32Hz','64Hz','128Hz'};
end
%% Plot TTF
% fit with 4th order polynomial
x = 1:1:length(means);
Pm = polyfit(x,means,4);
Ps = polyfit(x,SEMs,4);
xx = 1:0.01:length(means);
mm = polyval(Pm,xx);
ss = polyval(Ps,xx);
% Plot
figure;fill([(xx),fliplr(xx)],[(mm+ss),fliplr(mm-ss)],...
    [0.75 0.75 0.75],'LineWidth',1,'EdgeColor','none');
hold on;
plot(xx,mm,'k','LineWidth',2); hold on;
plot([xTick(1) xTick(end)],[0 0],'k'); hold on;
plot(x,means,'.r','MarkerSize',20);
%figure;errorbar(mm,ss);
title([hemi ' ' ROI ' ' condName ' ' func],'FontSize',30)
xlabel('Frequency','FontSize',20);
ylabel('Percent Signal Change','FontSize',20);
xlim(xlims);
ylim(ylims);
ax = gca;
set(ax,'XTick',xTick);
set(ax,'XTickLabel',xLabels);
axis square;