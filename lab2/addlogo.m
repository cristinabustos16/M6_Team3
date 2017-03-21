%%%%% Add logo to image.

clearvars
close all

% Destination image:
des = imread('./Data/castle_int/0016_s.png');

% Logo:
logo = imread('onlylogo_cocacola.png');
% logo = imread('onlylogo_ibm.png');
% logo = imread('onlylogo_gadis.png');

% Point 1 (upper left corner):
% pdes1 = [216, 401];
pdes1 = [548, 375];

% Point 2 (upper right corner):
% pdes2 = [254, 405];
pdes2 = [665, 362];

% Point 3 (lower left corner):
% pdes3 = [214, 426];
pdes3 = [549, 421];

% Point 41 (lower right corner):
% pdes4 = [253, 428];
pdes4 = [667, 414];

% Show destination points:
figure()
imshow(des)
hold on
plot(pdes1(1), pdes1(2), '*r')
plot(pdes2(1), pdes2(2), '*g')
plot(pdes3(1), pdes3(2), '*b')
plot(pdes4(1), pdes4(2), '*y')

% Source points (corners of the logo):
psrc1 = [0, 0];
psrc2 = [size(logo,2), 0];
psrc3 = [0, size(logo,1)];
psrc4 = [size(logo,2), size(logo,1)];

% Compute the homography between the source and destination points:
x1 = [psrc1', psrc2', psrc3', psrc4']; % cartesian
x1 = [x1; ones(1,4)]; % homogeneous
x2 = [pdes1', pdes2', pdes3', pdes4']; % cartesian
x2 = [x2; ones(1,4)]; % homogeneous
H = homography2d(x1, x2);

corners = [1, size(des,2), 1, size(des,1)];
logo_trans = apply_H_v2(logo, H, corners);

figure()
imshow(logo_trans)

% mask(:,:,1) = uint8(logo_trans(:,:,3) < 100); % ad-hoc para ibm
mask(:,:,1) = uint8(logo_trans(:,:,1) < 5); % ad-hoc para gadis
% mask(:,:,1) = uint8(logo_trans(:,:,1) < 100); % ad-hoc para cocacola
mask(:,:,2) = mask(:,:,1);
mask(:,:,3) = mask(:,:,1);
% mask = uint8(logo_trans == 0);

newimage = logo_trans + des .* mask;
figure()
imshow(newimage)
