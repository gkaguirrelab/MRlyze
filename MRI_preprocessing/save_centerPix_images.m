function save_centerPix_images(inputImg,outputImg,centerPix,dotRadius)

% Saves a red dot on the 'centerPix' in an image
%
%   Written by Andrew S Bock Dec 2015

%% set defaults
if ~exist('dotRadius','var')
    dotRadius = 20; % radius of 2 pixels
end
%% load image
Im = imread(inputImg);
dims = size(Im);
%% Find distance from center pixel, color red if within dotRadius
dists = nan(dims(1)*dims(2),1);
for i = 1:dims(1)
    for j = 1:dims(2)
        dst = sqrt( (i - centerPix(1)).^2 + (j - centerPix(2)).^2 );
        if dst <= dotRadius
            Im(i,j,:) = [255 0 0];
        end
    end
end
%% Save image
imwrite(Im,outputImg,'png');
