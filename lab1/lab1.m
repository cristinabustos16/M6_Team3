%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 1: Image rectification


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Applying image transformations

% ToDo: create the function  "apply_H" that gets as input a homography and
% an image and returns the image transformed by the homography.
% At some point you will need to interpolate the image values at some points,
% you may use the Matlab function "interp2" for that.


%% 1.1. Similarities
I=imread('Data/0005_s.png'); % we have to be in the proper folder

% ToDo: generate a matrix H which produces a similarity transformation

% Similarity transformation: let x and x??? be homogeneous coordinates. 
% x' = H_s x, where H_s = (sR, t; 0^T, 1): s is an isotropic scaling factor, 
% R is a rotation matrix (orthogonal matrix) and t is a translation vector.

scalingFactor = 1;
rotationAngle = 45 / 360 * 2 * pi;
translationX = 0;
translationY = 0;

H=[ scalingFactor * cos(rotationAngle)   scalingFactor * -sin(rotationAngle)    translationX;
    scalingFactor * sin(rotationAngle)   scalingFactor *  cos(rotationAngle)    translationY;
    0                                    0                                      1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));


%% 1.2. Affinities

% ToDo: generate a matrix H which produces an affine transformation

A = [ 1 0.3;
     -1 0.6];
translationX = 2;
translationY = 3;

H = [ A(1,1)   A(1,2)    translationX;
      A(2,1)   A(2,2)    translationY;
      0        0         1];
I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation

[U,D,V] = svd(A);
RTheta = U*V';  % U --> buen resultado visual
RThetaTransform = [RTheta(1,1) RTheta(1,2) 0;
                   RTheta(2,1) RTheta(2,2) 0; 
                   0           0           1];
RPhi = V';
RPhiTransform = [RPhi(1,1) RPhi(1,2) 0;
                 RPhi(2,1) RPhi(2,2) 0; 
                 0         0         1];
lambda1 = D(1,1);
lambda2 = D(2,2);
scaleTransform = [lambda1 0       0;
                  0       lambda2 0;
                  0       0       1];
translationTransform = [1 0 translationX;
                        0 1 translationY;
                        0 0 1];
 
% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above

HDecomposed = translationTransform * (RThetaTransform * RPhiTransform' * scaleTransform * RPhiTransform);
error_threshold = 1e-10;
if sum(sum(round(H * 100) / 100 - round(HDecomposed * 100) / 100)) == 0
    fprintf('yeah! the product of the four previous transformations produces the same matrix H as above');
else
    fprintf('ohhhh! the product of the four previous transformations does not produce the same matrix H as above');
end

% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before
 
I2Decomposed = apply_H(I, HDecomposed);
fig = figure(2);
subplot(1,3,1); imshow(I); title('Original image');
subplot(1,3,2); imshow(I2); title('Affine transformation with H');
subplot(1,3,3); imshow(I2Decomposed); title('Affine transformation with H decomposed');


%% 1.3 Projective transformations (homographies)

% ToDo: generate a matrix H which produces a projective transformation

A = [0.6 0.1;
     0.3 0.4];
translationX = 2;
translationY = 3;
v1 = 0.0001;
v2 = 0.0002;
v = 1;

H = [A(1,1)   A(1,2)    translationX;
     A(2,1)   A(2,2)    translationY;
     v1       v2        v];
I2 = apply_H(I, H);

v = -1;
H = [A(1,1)   A(1,2)    translationX;
     A(2,1)   A(2,2)    translationY;
     v1       v2        v];
I3 = apply_H(I, H);

fig = figure(3);
subplot(1,3,1); imshow(I); title('Original image');
subplot(1,3,2); imshow(I2); title('Projective transformation');
subplot(1,3,3); imshow(I3); title('Projective transformation mirrored');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification

% choose the image points
I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');
affineRectification(I, A, 424, 240, 712, 565);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.

% We will use the same image, points and lines as in the previous sections.
% l1 and l3 are perpendicular, the same as l2 and l4.

% We actually use the rectified lines:

B = [l1_a(1) * l3_a(1), l1_a(1) * l3_a(2) + l1_a(2) * l3_a(1), l1_a(2) * l3_a(2);
     l2_a(1) * l4_a(1), l2_a(1) * l4_a(2) + l2_a(2) * l4_a(1), l2_a(2) * l4_a(2)];
[~, ~, V] = svd(B);
s = V(:,end); % I hope this is the null vector of B.
S = [s(1), s(2);
     s(2), s(3)];

% Compute the upper triangular matrix using the Cholesky factorization:
A = chol(S, 'lower');
% A = A / sqrt(det(A));

% Matrix for the metric rectification:
H_a = [A(1,1), A(1,2), 0;
       A(2,1), A(2,2), 0;
       0     , 0     , 1];

H = H_p * H_a;

% I_s = apply_H(I, H);
% figure; imshow(uint8(I_s));
% 
% figure()
% subplot(1,3,1)
% imshow(uint8(I))
% subplot(1,3,2)
% I_a = apply_H(I, H_p);
% imshow(uint8(I_a))
% subplot(1,3,3)
% imshow(uint8(I_s))

% Angulos entre lineas:
l1_s = inv(H)' * l1;
l2_s = inv(H)' * l2;
l3_s = inv(H)' * l3;
l4_s = inv(H)' * l4;

l1_s = inv(H_a)' * l1_a;
l2_s = inv(H_a)' * l2_a;
l3_s = inv(H_a)' * l3_a;
l4_s = inv(H_a)' * l4_a;

compute_angle(l1_s, l3_s)
compute_angle(l2_s, l4_s)




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 4. OPTIONAL: Metric Rectification in a single step
% % Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)
% 

%% 5. OPTIONAL: Affine Rectification of the left facade of image 0000
I = imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');
affineRectification(I, A, 493, 186, 48, 508);

%% 6. OPTIONAL: Metric Rectification of the left facade of image 0000
% 

%% 7. OPTIONAL: Affine Rectification of the left facade of image 0001
I = imread('Data/0001_s.png');
A = load('Data/0001_s_info_lines.txt');
affineRectification(I, A, 614, 159, 645, 541);
 
%% 8. OPTIONAL: Metric Rectification of the left facade of image 0001
% 

