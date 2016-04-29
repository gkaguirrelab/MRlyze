function [varexp,params,sorted_templates] = find_best_template(templateType,tdir,hemi,psi,FCx,FCy,fitType)

% Finds the best pRF model template
%
%   Usage:
%   [varexp,params,sorted_templates] = find_best_template(templateType,tdir,hemi,psi,FCx,FCy)
%
%   Output:
%   template parameters, sorted by variance explained
%
%   Defaults (coarse):
%   psi = [-0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
%   FCx = [-0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6];
%   FCy = [-0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5];
%
%   Defaults (fine):
%   psi = [-0.066 -0.033 0.000 0.033 0.066];
%   FCx = [-0.066 -0.033 0.000 0.033 0.066];
%   FCy = [-0.066 -0.033 0.000 0.033 0.066];
%
%   Written by Andrew S Bock Sep 2015

%% Set defaults
% template parameters
if ~exist('psi','var') || isempty(psi)
    if strcmp(templateType,'coarse')
        %psi = [-0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
        psi = [0.4 0.5 0.6 0.7 0.8 0.9 1.0];
    elseif strcmp(templateType,'fine')
        psi = [-0.066 -0.033 0.000 0.033 0.066];
    end
end
if ~exist('FCx','var') || isempty(FCx)
    if strcmp(templateType,'coarse')
        %FCx = [-0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6];
        FCx = [-0.45 -0.35 -0.25 -0.15 -0.05 0.05 0.15];
    elseif strcmp(templateType,'fine')
        FCx = [-0.066 -0.033 0.000 0.033 0.066];
    end
    
end
if ~exist('FCy','var') || isempty(FCy)
    if strcmp(templateType,'coarse')
        %FCy = [-0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5];
        FCy = [-0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1];
    elseif strcmp(templateType,'fine')
        FCy = [-0.066 -0.033 0.000 0.033 0.066];
    end
end
if ~exist('fitType','var')
   fitType = 'V1'; % 'V1 = V1<->V2, V1<->V3; 'V2V3' =  V1<->V2, V1<->V3, AND V2<->V3
end
%% Pull out the variance explained
varFiles = listdir(fullfile(tdir,[hemi '*varexp.txt']),'files');
% remove the 'pRF' varexp files, if they exist
pRFfiles = strfind(varFiles,'pRF');
pRFind = [];
for i = 1:length(pRFfiles);
    if ~isempty(pRFfiles{i})
        pRFind = [pRFind i];
    end
end
varFiles(pRFind) = [];
% remove the 'anat' varexp files, if they exist
anatfiles = strfind(varFiles,'anat');
anatind = [];
for i = 1:length(anatfiles);
    if ~isempty(anatfiles{i})
        anatind = [anatind i];
    end
end
varFiles(anatind) = [];
% load in the variance explained text files
if strcmp(fitType,'V1')
    vals = zeros(length(varFiles),4);
elseif strcmp(fitType,'V2V3')
    vals = zeros(length(varFiles),6);
end
for i = 1:length(varFiles);
    vals(i,:) = load(fullfile(tdir,varFiles{i}));
end
%% Extract the variance explained
if strcmp(fitType,'V1')
    subvals = vals(:,1:2); % V1-V2, V1-V3
elseif strcmp(fitType,'V2V3')
    subvals = vals(:,1:3); % V1-V2, V1-V3, V2-V3
end
sumvals = nansum(subvals,2);
%% Sort the values, output the template parameters
[varexp,sortind] = sort(sumvals);
varexp = flipud(varexp);
sortind = flipud(sortind);
for i = 1:length(sortind)
    dotinds = strfind(varFiles{sortind(i)},'.');
    params(i).psi = psi(str2double(varFiles{sortind(i)}(dotinds(1)+1:dotinds(2)-1)));
    params(i).FCx = FCx(str2double(varFiles{sortind(i)}(dotinds(2)+1:dotinds(3)-1)));
    params(i).FCy = FCy(str2double(varFiles{sortind(i)}(dotinds(3)+1:dotinds(4)-1)));
    sorted_templates{i} = varFiles{sortind(i)};
end