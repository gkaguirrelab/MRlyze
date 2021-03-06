function firFig = FIR_plot(means, SEMs, ROI, condition, hemi, func, subj_name, runNums, xlims,ylims,xTick,xLabels)
% firFig = FIR_plot(means, SEMs, ROI, condition, hemi, func, subj_name, runNums, ...
% xlims,ylims,xTick,xLabels)
%
% Plots FIR.
%
% 
%
% Input arguments:
% ================
%
%   means : means values to plot;
%   SEMs : SEMs value to plot;
%   ROI : Region of Interest (for plot title);
%   condition: condition name (for plot title);
%   hemi : hemisphere (for plot title);
%   func : functional volume (for plot title);
%   subj_name : subject name (for plot title);
%   runNums : number of valid runs;
%   xlims : limit on the X axis;
%   ylims : limit on the Y axis;
%   xTick : ticks one X axix;
%   xLabels : labels on the X axis.
% 
% 
% Usage:
% ======
%   this function is used within the processing function
%   'MelanopsinMR/analysis/MelanopsinMR_PostFeatAnalysis.m'
% 
%
%
% 3/1/2016    gf, ms Based on plot_TTF, Written by Andrew S Bock Sep 2015
% 7/2/2016    ms Commented.
% 7/6/2016    gf updated y axis limit to avoid off-chart plots

%% set defaults
if ~exist('condName','var')
    CondName = '';
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
if ~exist('subj_name','var')
    subj_name = '';
end
if ~exist('runNums','var')
    runLabel = false;
else
    runLabel = true;
end
if ~exist('xlims','var')
    xlims = [-1 15];
end
if ~exist('ylims','var')
    ylims = [-0.5 1.4];
end
if ~exist('xTick','var')
    xTick = [0 1 2 3 4 5 6 7 8 9 10 11 12 13];
end
if ~exist('xLabels','var')
    xLabels = xTick;
end

%% Plot FIR with errorbars
x = 0:1:(length(means)-1);

firFig = figure('units','normalized','position',[0 0 1 1]);
plot([xlims(1) xlims(end)],[0 0],'k'); hold on;
h = plot(x, means, ':k'); hold on;
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
g = errorbar(x,means,SEMs, '.k','MarkerSize',16); hold on;
set(get(get(g,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
p = plot(x, means, '.r', 'MarkerSize',16);
legend (p, ['Mean\pmSEM']);
title({subj_name, [condition ' - ' hemi ' ' ROI ' ' func]},'Interpreter','none')
xlabel('Time [sec]');
ylabel('Amplitude [% signal change]');
xlim(xlims); ylim(ylims);
ax = gca;

if runLabel
str = ['Valid runs = ',num2str(length(runNums))];
text(0, -0.4, str)
end

set(ax,'XTick',xTick);
set(ax,'XTickLabel',xLabels);