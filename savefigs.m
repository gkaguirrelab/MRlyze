function savefigs(savetype,filename)
% This function allows you to quickly save all currently open figures with
% a custom filename for each in multiple formats.  To use the function
% simply call savefigs with no arguments, then follow the prompts
%
% Upon execution this function will one-by-one bring each currently open
% figure to the foreground.  Then it will supply a text prompt in the main
% console window asking you for a filename.  It will save that figure to
% that filename in the .fig (using saveas) and .tif (using print) formats.
%
% The formats that it saves in can be changed by commenting out or adding
% lines below.
%
% Copyright 2010 Matthew Guidry
% matt.guidry ATT gmail DOTT com  (Email reformatted for anti-spam)
% Updated 2012 Andrew Bock
% using 'saveas' for figures resulted in poor resolution/dimensions
% now uses 'print'

if ~exist('savetype','var')
    disp('You must specify a file type to save, e.g. ''pdf'',''png'', or ''eps''')
    return
end
hfigs = get(0, 'children') ;                         %Get list of figures

for m = 1:length(hfigs)
    figure(hfigs(m))                                %Bring Figure to foreground
    %     screen_size = get(0, 'ScreenSize');                             %Make Figure fullscreen
    %     full = figure(hfigs(m));                                        %Make Figure fullscreen
    %     set(full, 'Units','normalized','Position',[.5 1 .5 1]); %Make Figure fullscreen
    if ~exist('filename','var')
        filename = input('Filename? (0 to skip)\n', 's');%Prompt user
    end
    if strcmp(filename, '0')                        %Skip figure when user types 0
        continue
    else
        set(hfigs(m),'PaperPositionMode','auto')
        %set(gcf, 'PaperUnits', 'normalized');
        %x_width=10 ;y_width=10;
        %set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
        %print(filename,'-dtiff','-r300') % save at 300 dpi
        %saveas(hfigs(m), [filename '.fig']) %Matlab .FIG file
        %set(gcf, 'PaperPosition', [0 0 20 12]); %Position plot at left hand corner with width 20 and height 12.
        if strcmp(savetype,'png')
            saveas(hfigs(m),filename,'png') %Save as png
        elseif strcmp(savetype,'eps');
            saveas(hfigs(m),filename,'eps') %Save as eps
        elseif strcmp(savetype,'pdf');
            set(hfigs(m),'Units','Inches');
            pos = get(hfigs(m),'Position');
            set(hfigs(m),'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
            saveas(hfigs(m),filename,'pdf') %Save as pdf
        elseif strcmp(savetype,'all')
            saveas(hfigs(m),filename,'png') %Save as png
            saveas(hfigs(m),filename,'eps') %Save as eps
            set(hfigs(m),'Units','Inches');
            pos = get(hfigs(m),'Position');
            set(hfigs(m),'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
            saveas(hfigs(m),filename,'pdf') %Save as pdf
        end
    end
    clear filename
end