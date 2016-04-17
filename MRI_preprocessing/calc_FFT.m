function [out] = calc_FFT(in,fs,makePlot)

% calculates the frequency spectrum on an input signal
%
%   in = input signal
%   fs = sampling frequency (e.g. 50 Hz)
%
%   Written by Andrew S Bock Dec 2015

%% Create the output frequency spectrum
out = abs(fftshift(fft(detrend(in,'constant'))));

%% Plot the frequencies
if makePlot
    figure;
    subplot(1,2,1);
    plot(in);hold on
    axis tight;
    df = fs/length(in); % calculate the delta
    x = -fs/2:df:(fs/2 - df); % make the x axis values
    subplot(1,2,2);
    plot(x,out);
    axis tight;
end

