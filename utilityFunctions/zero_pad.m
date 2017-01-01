function zero_pad(indir)

% Zero pads .mat file names such that each mat file name has 4 digits
%   Originally written to convert files in /STIMMATRIX folder in the
%   subject directories of the /jet/aguirre/MF_Retinotopy/LongRetinotopy
%   project directory
%
%   Written by Andrew S Bock 2015

mats = listdir(fullfile(indir,'*.mat'),'files');
for mm = 1:length(mats);
    infile = mats{mm};
    TRstr = strfind(infile,'TR');
    dotmatstr = strfind(infile,'.mat');
    innum = infile(TRstr+2:dotmatstr-1);
    outnum = sprintf('%04d',str2double(innum));
    outfile = [infile(1:TRstr+1) outnum infile(dotmatstr:end)];
    movefile(infine,outfile);    
end