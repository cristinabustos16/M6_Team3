%%%% Create a mask for the CVC logo.

clearvars
close all

logo = imread('logo_cvc.png');

mask = uint8(logo(:,:,2) > 180);

onlylogo(:,:,1) = logo(:,:,1) .* (1 - mask);
onlylogo(:,:,2) = logo(:,:,2) .* (1 - mask);
onlylogo(:,:,3) = logo(:,:,3) .* (1 - mask);

mask_open = imopen(mask, strel('disk',2));
onlylogo_open(:,:,1) = logo(:,:,1) .* (1 - mask_open);
onlylogo_open(:,:,2) = logo(:,:,2) .* (1 - mask_open);
onlylogo_open(:,:,3) = logo(:,:,3) .* (1 - mask_open);

figure()
subplot(2,2,1)
imshow(logo)
subplot(2,2,2)
imshow(mask, [0, 1])
subplot(2,2,3)
imshow(onlylogo)
subplot(2,2,4)
imshow(onlylogo_open)


% Save the onlylogo (I like it more than with the opening) and the mask:
imwrite(onlylogo, 'onlylogo_cvc.png');
mask = mask * 255;
imwrite(mask, 'mask_logo_cvc.png');

