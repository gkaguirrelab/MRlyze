function filterMotionMatFiles(parFile,MatDir,outDir,TR)

% Filters motion parameters at twice the Nyquist frequency, saves as MAT 
%   files as expected by FSL's mclirt
%
%   Usage:
%   filterMotionMatFiles(parFile,MatDir,outDir,TR)
%
%   note:
%   FSL is strange for these params:
%   https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1112&L=FSL&D=0&P=7568
%
%   Thus, we use the MAT files for the translations, and the .par file for 
%   the rotations.
%
%   note:
%   TR is in seconds
%
%   Based on scripts by Marta Vidorreta Díaz de Cerio
%   Written by Andrew S Bock Apr 2016

%% Read in motion params
% parfile for rotations
fid = fopen(parFile,'r');
Amotion = (fscanf(fid, '%g %g %g %g %g %g\r\n',[6 inf]))';
fclose(fid);
mc_rotations = Amotion(:,1:3);
% MAT files for translations
matFiles = listdir(MatDir,'files');
n_scans = length(matFiles);
mc_translations = zeros(n_scans,3);
for i = 1:n_scans
    fid=fopen(fullfile(MatDir,matFiles{i}),'r');
    MAT = fscanf(fid,'%g %g %g %g\r\n',[4 inf])';
    fclose(fid);
    mc_translations(i,1:3) = (MAT(1:3,4))';
end
motParams = [mc_rotations mc_translations];
%% Filter the Nyquist frequency
inSignal = motParams;
filtType = 'low';
sampT = TR;
cutoffHzlow = -inf;
cutoffHzhigh = 1/(4*TR); % Twice the Nyquist frequency
[filtMotParams] = filter_signal(inSignal,filtType,sampT,cutoffHzlow,cutoffHzhigh);
%% Loop through motion params
ct = 0;
for j = 1:size(filtMotParams,1);
    ct = ct + 1;
    % Create matrix for this TR
    [outMat] = convertMotion2Mat(filtMotParams(j,:));
    % Save the MAT file
    fid=fopen(fullfile(outDir,['MAT_' sprintf('%04d',ct)]),'w');
    fprintf(fid,'%1.6f %1.6f %1.6f %1.6f\r\n',outMat);
    fclose(fid);
end