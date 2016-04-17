function convert_image2surf(inputImg,outputNifti,dataEcc,dataPol,datatAreas,centerPix,display)

% Converts an input image into a surface overlay
%
%   Usage:
%   convert_image2surf(inputImg,outputNifti,dataEcc,dataPol,centerPix,display)
%
%   Inputs:
%   inputImage = '/some/image/file.png';
%   outputNifti = '/some/nifti/file.nii.gz';
%   dataEcc = '/some/Ecc/file/or/data/matrix';
%   dataPol = '/some/Pol/file/or/data/matrix'; % note, must be in RADIANS!
%   dataAreas = '/some/Areas/file/or/data/matrix';
%   centerPix = [600 960]; % default - assumes 1920 x 1200 resolution
%   (image loads rows x column)
%   display.distance = 180.2; % default centimeters - SC7T UPenn
%   display.resolution = [1920 1080]; % default screen [width height] - SC7T UPenn
%   display.width = 69.7347; % default centimeters - SC7T UPenn
%   display.heigth = 39.2257; % default centimeters - SC7T UPenn
%
%   Written by Andrew S Bock Dec 2015

%% set defaults
if ~exist('display','var')
    display.distance = 180.2; % default centimeters - SC7T UPenn
    display.resolution = [1920 1080]; % screen [width height] - SC7T UPenn
    display.width = 69.7347; % default centimeters - SC7T UPenn
    display.heigth = 39.2257; % default centimeters - SC7T UPenn
end
%% Load input image
Im = imread(inputImg);
if ~exist('centerPix','var')
    centerPix = round([size(Im,1)/2 size(Im,2)/2]); % default - assumes 1920 x 1200 resolution
end
% Set to pixel coordinates (using fixation 'centerPix')
Imx = (1:size(Im,1)) - centerPix(1);
Imy = (1:size(Im,2)) - centerPix(2);
%% Load Eccentricity
if ischar(dataEcc)
    if exist(dataEcc,'file')
        [~,~,ext] = fileparts(dataEcc);
        if strcmp(ext,'.gz')
            tmp = load_nifti(dataEcc);
            srfEcc = tmp.vol;
        elseif strcmp(ext,'.mgh') || strcmp(ext,'.mgz')
            srfEcc = load_mgh(dataEcc);
        end
    else
        error('data input does not correspond to an existing file');
    end
else
    srfEcc = double(dataEcc);
end
%% Load Polar Angle
if ischar(dataPol)
    if exist(dataPol,'file')
        [~,~,ext] = fileparts(dataPol);
        if strcmp(ext,'.gz')
            tmp = load_nifti(dataPol);
            srfPol = tmp.vol;
        elseif strcmp(ext,'.mgh') || strcmp(ext,'.mgz')
            srfPol = load_mgh(dataPol);
        end
    else
        error('data input does not correspond to an existing file');
    end
else
    srfPol = double(dataPol);
end
%% Load Areas
if ischar(datatAreas)
    if exist(datatAreas,'file')
        [~,~,ext] = fileparts(datatAreas);
        if strcmp(ext,'.gz')
            tmp = load_nifti(datatAreas);
            srfAreas = tmp.vol;
        elseif strcmp(ext,'.mgh') || strcmp(ext,'.mgz')
            srfAreas = load_mgh(datatAreas);
        end
    else
        error('data input does not correspond to an existing file');
    end
else
    srfAreas = double(datatAreas);
end
%% Convert surface ecc and angle values to pixel space
srfEccPix = angle2pix(display,srfEcc);
[tmpy,tmpx] = pol2cart(srfPol,srfEccPix);
srfx = round(tmpx);
srfy = round(tmpy);
%% Replace surface values with image values
newsurf = nan(length(srfx),3);
V1V3ind = find(abs(srfAreas)<=3);
newsurf(V1V3ind,:) = repmat([128 128 128],length(V1V3ind),1); % set V1-V3 to gray
for i = 1:length(srfx)
    ImxInd = find(Imx==srfx(i),1);
    ImyInd = find(Imy==srfy(i),1);
    if ~isempty(ImxInd) && ~isempty(ImyInd);
        newsurf(i,:) = Im(ImxInd,ImyInd,:);
    end
end
%% Save output file
tmp.vol = newsurf./255;
tmp.vol = reshape(tmp.vol,size(tmp.vol,1),1,1,3);
tmp.dim(5) = 3;
save_nifti(tmp,outputNifti);