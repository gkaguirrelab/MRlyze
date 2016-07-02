function FIR_multiplot (output_dir,dataFiles, dataExt, legendTexts, titleText, saveName)
% FIR_multiplot (output_dir,dataFiles, dataExt, legendTexts, titleText, saveName)
%
%  plots multiple FIR responses on the same graph. Gets data values from either
%  .fig files or .csv files.
%
% Example:
%
% output_dir ='/Users/giulia/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/MelanopsinMR/MelPulses_400%/All_subjects';
% dataFiles = { ...
%     '/data/jag/MELA/HERO_asb1/032416/FIR_figures/mh_V1_MaxMelPulse_FIR_raw.fig' ...
%     '/data/jag/MELA/HERO_aso1/032516/FIR_figures/mh_V1_MaxMelPulse_FIR_raw.fig' ...
%     '/data/jag/MELA/HERO_gka1/033116/FIR_figures/mh_V1_MaxMelPulse_FIR_raw.fig' ...
%     '/data/jag/MELA/HERO_mxs1/040616/FIR_figures/HERO_mxs1_MaxMelPulse_mh_V1_FIR_raw.fig' ...
%     };
% dataExt = 'fig' ;
% legendTexts = { ...
%     'HERO_asb1' ...
%     'HERO_aso1' ...
%     'HERO_gka1' ...
%     'HERO_mxs1' ...
%     };
% titleText = 'All subjects MaxMel_mh_V1';
% saveName = 'HEROES_MaxMel_mh_V1' ;
%


%% set axes
figure('units','normalized','position',[0 0 1 1]);
hold on
xlabel('Time [sec]');
ylabel('Signal change [%]');
title (titleText,'Interpreter','none')
xlims = [-1 15];
ylims = [-0.5 1.2];
xTick = [0 1 2 3 4 5 6 7 8 9 10 11 12 13];
xLabels = xTick;
h = plot([xlims(1) xlims(end)],[0 0],'k');
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
xlim(xlims); ylim(ylims);
ax = gca;
set(ax,'XTick',xTick);
set(ax,'XTickLabel',xLabels);

for ff = 1:length(dataFiles)
    switch dataExt
        case 'fig' % Get datapoints from fig files
            H = open (dataFiles{ff});
            D=get(gca,'Children');
            YData=get(D,'YData');
            y(ff) = YData(1);
            close (H);
            
            dataP = y(ff);
            dataP = transpose (dataP{:});
        case 'csv' % Get datapoints from csv files
            y(:,ff) = csvread(dataFiles{ff});
            dataP = y(:,ff);
    end
    dataL = length(dataP);
    x = 0:1:dataL-1;
    plot(x,dataP,'o-');
    axis square;
    legendInfo{ff} = (legendTexts{ff});
end
legend (legendInfo, 'Interpreter','none');

% Save pdf in output dir
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

adjustPlot(gcf);
saveas(gcf, fullfile(output_dir, saveName), 'pdf');
close all;