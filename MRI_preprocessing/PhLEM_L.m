function [x,events] = PhLEM_L(TS,Hz,HRrate)
% Filter and smoothing function for LOW frequencies
if strcmp(HRrate,'low')
    delta = Hz*1.5;
    Wn = Hz*0.3;
else
    delta = Hz*1.5; %Default is Hz*1.5
    Wn = Hz*0.35; %Default is Hz*0.4
end
peak_rise = 0.5;
% Transform via ABS
x = abs(TS);
x = x - mean(x);
% Use Gaussian Filter
x = smooth_kernel(x(:),Wn);
x = x';
% -- Schlerf (Jul 24, 2008)
% Take a first pass at peak detection:
[mxtb,mntb] = peakdet(x,1e-14);
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
[mxtb,mntb] = peakdet(x,peak_rise*peak_resp);
events = zeros(size(TS));
peaks = mxtb(:,1);
dpeaks = diff(peaks);
kppeaks = find(dpeaks > delta);
newpeaks = peaks([1 kppeaks'+1]);
events(newpeaks) = 1;
% DEBUG CATCH:
if length(newpeaks) < 5;
    delta = Hz*1.275;
    [mxtb,mntb] = peakdet(x,1e-14);
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
    [mxtb,mntb] = peakdet(x,peak_rise*peak_resp);
    events = zeros(size(TS));
    peaks = mxtb(:,1);
    dpeaks = diff(peaks);
    kppeaks = find(dpeaks > delta);
    newpeaks = peaks([1 kppeaks'+1]);
    events(newpeaks) = 1;
    disp('Low Freq peaks hard to detect, switching delta')
    if length(newpeaks) < 5;
        warning('Program Crash: Low Freq peaks still hard to detect; Bad Data?')
    end
    %   keyboard;
end