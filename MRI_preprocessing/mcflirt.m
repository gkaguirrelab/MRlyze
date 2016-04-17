function mcflirt(in,out)
%   Calls FSL's "mcflirt" command to motion correct epi data. This function
%   will call "mcflirt -in <in> -o <out> -stats -mats -plots -report
%
%   Usage;
%   mcflirt(in,out)
%
%   Written by Andrew S Bock Sept 2014

if ~exist('in','var')
    error('no "in" image defined');% must define an input image
end
if ~exist('out','var')
    error('no "out" image defined');% must define a output image
end
[~,inname,~] = fileparts(in);
[~,outname,~] = fileparts(out);
%%
system(['mcflirt -in ' inname ' -o ' outname ' -stats -mats -plots -report']);