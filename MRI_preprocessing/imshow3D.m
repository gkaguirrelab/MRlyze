function img2dshow = imshow3D(input_vol)
% Converts 3D image into a 2D mosaic
%   Displays with imshow with no blank spaces between image slices.
%
% Usage:
%   img2dshow = imshow3D(input_vol)
%
% Originally written by Marta Vidorreta Diaz de Cerio
%   Updated by Andrew S Bock Feb 2015

%% Load in 3D nifti
disp(['Loading ' input_vol '...']);
tmp = load_nifti(input_vol);
input_img = tmp.vol;
disp('done.');
%% Create 2D mosaic
disp('Creating 2D mosaic...');
img = permute(input_img,[2 1 3]);
n_slices = size(img,3);
n_cols = ceil(sqrt(n_slices));
n_rows = ceil(n_slices/n_cols);
if n_rows*n_cols > n_slices
    n_dif = n_rows*n_cols - n_slices;
    img(:,:,end+1:end+n_dif) = 0;
end
img2d = flipud(reshape(img, size(img,1), size(img,2)*size(img,3)));
img2dshow = zeros(size(img,1)*n_rows, size(img,2)*n_cols);
for irow = 1:n_rows
    I1r = 1 + (irow-1)*size(img,1);
    I2r = I1r + size(img,1) - 1;
    I1c = 1 + (irow-1)*n_cols*size(img,2);
    I2c = I1c + n_cols*size(img,2) - 1;
    img2dshow(I1r:I2r, :) = img2d(:, I1c:I2c);
end
disp('done.');
%% Plot the mosaic
disp(['Plotting 2D mosaic for ' input_vol]);
figure;
imshow(img2dshow, [min(img2dshow(:)) max(img2dshow(:))]);
colorbar;