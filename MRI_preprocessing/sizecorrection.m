function [outputsignal] = sizecorrection(signal,newlength)

x = 1:length(signal);
shifter = length(signal)/newlength;

outputsignal = zeros (1,newlength);

currentspot = 1;
for cc=1:newlength
    
    if currentspot > length(signal)
        currentspot = length(signal);
    else
        outputsignal(1,cc) = interp1(x,signal,currentspot,'linear');
        currentspot = currentspot + shifter;
    end
end