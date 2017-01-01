function outImage = convert_surf2image(inVol,inFrames,eccVol,polVol,axLim,vArea,matSize,sigScale)

% Converts input surface to visual field space
%
%   Usage:
%   outImage = convert_surf2image(inVol,inFrames,eccVol,polVol,roiInd,axLim,vArea)
%
%   Written by Andrew S Bock Jun 2016

%% set defaults
if ~exist('axLim','var')
    axLim = 6.2;
end
if ~exist('vArea','var')
    vArea = 'V1';
end
if ~exist('matSize','var')
    matSize = 101;
end
%% Threshold based on axis limit
badind = eccVol>axLim;
eccVol(badind) = [];
polVol(badind) = [];
%% Convert polar to cartesian
[x,y] = pol2cart(polVol,eccVol);
y = -y; % flip y-axis

%% Create image mesh
centerMat = round(matSize/2);
[meshX,meshY] = meshgrid((1:matSize)-centerMat,(1:matSize)-centerMat);
matDist = sqrt( meshX.^2 + meshY.^2 );
polMat = cart2pol(meshX,meshY);
%% Log scale the image
eccMat = 10.^(matDist * log10(axLim)/(centerMat-1)); % log scale ecc
[matX,matY] = pol2cart(polMat,eccMat); % X and Y, log scale
%% Create sigma size based on visual area
Sig = rf_ecc(eccMat,vArea);
Sig = Sig * sigScale;
repSig = single(repmat(Sig,1,1,length(x)));
%% Create Gaussian weights for each point in the mesh
disp('Creating Gaussian weights...');
repmatX = single(repmat(matX,1,1,length(x)));
repmatY = single(repmat(matY,1,1,length(x)));
repx = single(repmat(x,1,matSize,matSize));
repx = permute(repx,[2,3,1]);
repy = single(repmat(y,1,matSize,matSize));
repy = permute(repy,[2,3,1]);
allDist = sqrt( (repmatY - repy).^2 + (repmatX - repx).^2);
clear repmatY repy repmatX repX
allGauss = exp(-(allDist.^2)./(2*repSig.^2));
clear allDist repSig 
flatGauss = single(reshape(allGauss,size(allGauss,1)*size(allGauss,2),size(allGauss,3)));
disp('done.');
%% Create out image(s)
disp('Creating output image...');
outImage = nan(matSize,matSize,length(inFrames));
ct = 0;
for f = inFrames
    ct = ct + 1;
    thisin = inVol(:,f);
    thisin(badind) = [];
    tmpOut = sum((flatGauss * thisin),2) ./ sum(flatGauss,2);
    tmpOut = reshape(tmpOut,[matSize matSize]);
    tmpOut(eccMat>axLim) = nan;
    outImage(:,:,ct) = tmpOut;
end
disp('done.');