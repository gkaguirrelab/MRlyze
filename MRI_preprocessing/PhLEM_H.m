function [x,events] = PhLEM_H(TS,Hz, HRrate)
% Filter and smoothing function for HIGH frequencies
if strcmp(HRrate,'low')
    delta = Hz/3.5;
    Wn = [1 10]/Hz;
else
    delta = Hz/6;
    Wn = [1 14]/Hz; % Default [ 1 12]; a range from 9 to 12 works for second number
end
peak_rise = 0.1;
%Transform via ABS
x = abs(TS);
x = x - mean(x);
% Use butter filter
if ~isempty(Wn)
    [b a] = butter(1,Wn);
    xbutter = filter(b,a,x);
end
% -- Schlerf (Jul 24, 2008)
% Take a first pass at peak detection:
mxtb  = peakdet(x,1e-14);
% set the threshold based on the 20th highest rather than the highest:
sorted_peaks = sort(mxtb(:,2));
if length(sorted_peaks) > 20
    peak_resp=sorted_peaks(end-20);
else
    if length(sorted_peaks) > 15
        peak_resp=sorted_peaks(end-15);
    else
        if length(sorted_peaks) > 10
            peak_resp=sorted_peaks(end-10);
        else
            if length(sorted_peaks) > 5
                peak_resp=sorted_peaks(end-5);
            else
                peak_resp = sorted_peaks(1);
            end
        end
    end
end
% And now do a second pass, more robustly filtered, to get the actual peaks:
mxtb  = peakdet(x,peak_rise*peak_resp);
events = zeros(size(TS));
peaks = mxtb(:,1);
dpeaks = diff(peaks);
kppeaks = find(dpeaks > delta);
newpeaks = peaks([1 kppeaks'+1]);
events(newpeaks) = 1;
% DEBUG CATCH:
if length(newpeaks) < 5;
    warning('Program Crash: High Freq peaks hard to detect; try different HRrate field?')
    keyboard;
end