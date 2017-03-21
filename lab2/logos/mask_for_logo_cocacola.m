%%%% Create a mask for the CVC logo.

clearvars
close all

logo = imread('logo_cocacola.png');

mask = uint8(logo(:,:,2) < 100);

onlylogo(:,:,1) = logo(:,:,1) .* mask;
onlylogo(:,:,2) = logo(:,:,2) .* mask;
onlylogo(:,:,3) = logo(:,:,3) .* mask;

figure()
subplot(2,2,1)
imshow(logo)
subplot(2,2,2)
imshow(mask, [0, 1])
subplot(2,2,3)
imshow(onlylogo)


% Save the onlylogo (I like it more than with the opening) and the mask:
imwrite(onlylogo, 'onlylogo_cocacola.png');
mask = mask * 255;
imwrite(mask, 'mask_logo_cocacola.png');

