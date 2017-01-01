function makeFEATDesignFile(templateFile,wildcards,parameters,outputFile)

%   Create a design file for running FSL FEAT, based on a template
%
%   Usage:
%   makeFEATDesignFile(templateFile,wildcards,parameters,outputFile)
%
%   Inputs
%   templateFile: path to .fsf file to be used as template
%       e.g. templateFile = ~/awesomeproject/FEAT/template.fsf
%   wildcards: cell array with one string in each row to be replaced by the
%       parameters values
%       e.g. wildcards = {'###SUBJID###';'###INPUTFILE###'}
%   parameters: cell array with each parameter coded as a string,
%       corresponding line-by-line to the wildcard cell array
%       e.g. parameters = {'sub001';'~/awesomeproject/sub001/BOLD1/f.nii.gz'}
%   outputFile: path to output .fsf file
%       e.g. outputFile = ~/awesomeproject/sub001/BOLD1/BOLD1.fsf
%
%   Written by Marcelo G Mattar April 2015
%

% Open template file
fin = fopen(templateFile,'rt');

% Open output file for edits
fout = fopen(outputFile,'wt');

% Replace wildcards with parameters
while(~feof(fin))
    s = fgetl(fin);
    for p = 1:length(wildcards)
        s = strrep(s,wildcards{p},parameters{p});
    end
    fprintf(fout,'%s\n',s);
    disp(s)
end
fclose(fin);
fclose(fout);