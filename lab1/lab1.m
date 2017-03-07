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

% % ToDo: generate a matrix H which produces a similarity transformation

% % Similarity transformation: let x and x??? be homogeneous coordinates. 
% % x' = H_s x, where H_s = (sR, t; 0^T, 1): s is an isotropic scaling factor, 
% % R is a rotation matrix (orthogonal matrix) and t is a translation vector.

% scalingFactor = 1;
% rotationAngle = 45 / 360 * 2 * pi;
% translationX = 0;
% translationY = 0;

% H=[ scalingFactor * cos(rotationAngle)   scalingFactor * -sin(rotationAngle)    translationX;
%     scalingFactor * sin(rotationAngle)   scalingFactor *  cos(rotationAngle)    translationY;
%     0                                    0                                      1];

% I2 = apply_H(I, H);
% figure; imshow(I); figure; imshow(uint8(I2));


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
if sum(sum(abs(H-HDecomposed))) < error_threshold
    fprintf('yeah! the product of the four previous transformations produces the same matrix H as above');
else
    fprintf('ohhhh! the product of the four previous transformations does not produce the same matrix H as above');
end

% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before
 
I2Decomposed = apply_H(I, HDecomposed);
fig = figure(2);
subplot(1,3,1); imshow(I); title('original');
subplot(1,3,2); imshow(I2); title('transformed with H');
subplot(1,3,3); imshow(I2Decomposed); title('transform with H decomposed');


% %% 1.3 Projective transformations (homographies)
% 
% % ToDo: generate a matrix H which produces a projective transformation
% 
% I2 = apply_H(I, H);
% figure; imshow(I); figure; imshow(uint8(I2));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 2. Affine Rectification


% % choose the image points
% I=imread('Data/0000_s.png');
% A = load('Data/0000_s_info_lines.txt');

% % indices of lines
% i = 424;
% p1 = [A(i,1) A(i,2) 1]';
% p2 = [A(i,3) A(i,4) 1]';
% i = 240;
% p3 = [A(i,1) A(i,2) 1]';
% p4 = [A(i,3) A(i,4) 1]';
% i = 712;
% p5 = [A(i,1) A(i,2) 1]';
% p6 = [A(i,3) A(i,4) 1]';
% i = 565;
% p7 = [A(i,1) A(i,2) 1]';
% p8 = [A(i,3) A(i,4) 1]';

% % ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
% l1 = cross(p1,p2);
% l2 = cross(p3,p4);
% l3 = cross(p5,p6);
% l4 = cross(p7,p8);

% % show the chosen lines in the image
% figure;imshow(I);
% hold on;
% t=1:0.1:1000;
% plot(t, -(l1(1)*t + l1(3)) / l1(2), 'r');
% plot(t, -(l2(1)*t + l2(3)) / l2(2), 'g');
% plot(t, -(l3(1)*t + l3(3)) / l3(2), 'b');
% plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% % % ToDo: compute the homography that affinely rectifies the image

% % compute the vanishing points, lecture 2 slide 16
% v1 = cross(l1, l2);
% v2 = cross(l3, l4);

% % compute the line at infinity
% l = cross(v1, v2);
% % normalize line at infinity
% l = l / max (l);

% % compute H based on the line at infinity, lecture2 slide 12
% H = [1 0 0; 0 1 0; l'];

% I2 = apply_H(I, H);
% figure; imshow(uint8(I2));

% % ToDo: compute the transformed lines lr1, lr2, lr3, lr4

% % Transpose of the inverse of H:
% H_inv_tras = (inv(H))';

% lr1 = H_inv_tras * l1;
% lr2 = H_inv_tras * l2;
% lr3 = H_inv_tras * l3;
% lr4 = H_inv_tras * l4;

% % show the transformed lines in the transformed image
% figure;imshow(uint8(I2));
% hold on;
% t=1:0.1:1000;
% plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'r');
% plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'g');
% plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'b');
% plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

% % ToDo: to evaluate the results, compute the angle between the different pair 
% % of lines before and after the image transformation

% % Lines before transformation in cartesian coordinates:
% l1c = [l1(1) / l1(3), l1(2) / l1(3)]';
% l2c = [l2(1) / l2(3), l2(2) / l2(3)]';
% l3c = [l3(1) / l3(3), l3(2) / l3(3)]';
% l4c = [l4(1) / l4(3), l4(2) / l4(3)]';

% % Angles of lines before transformation:
% ang_l1_l2_before = acos((l1c'*l2c)/sqrt((l1c'*l1c)*(l2c'*l2c)));
% ang_l3_l4_before = acos((l3c'*l4c)/sqrt((l3c'*l3c)*(l4c'*l4c)));
% fprintf('Angle of l1 and l2 before transformation = %f ?\n', ang_l1_l2_before)
% fprintf('Angle of l3 and l4 before transformation = %f ?\n', ang_l3_l4_before)

% % Lines after transformation in cartesian coordinates:
% l1rc = [lr1(1) / lr1(3), lr1(2) / lr1(3)]';
% l2rc = [lr2(1) / lr2(3), lr2(2) / lr2(3)]';
% l3rc = [lr3(1) / lr3(3), lr3(2) / lr3(3)]';
% l4rc = [lr4(1) / lr4(3), lr4(2) / lr4(3)]';

% % Angles of lines after transformation:
% ang_l1_l2_after = acos((lr1c'*lr2c)/sqrt((lr1c'*lr1c)*(lr2c'*lr2c)));
% ang_l3_l4_after = acos((lr3c'*lr4c)/sqrt((lr3c'*lr3c)*(lr4c'*lr4c)));
% fprintf('Angle of l1 and l2 after transformation = %f ?\n', ang_l1_l2_after)
% fprintf('Angle of l3 and l4 after transformation = %f ?\n', ang_l3_l4_after)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 3. Metric Rectification
% 
% %% 3.1 Metric rectification after the affine rectification (stratified solution)
% 
% % ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
% %       As evaluation method you can display the images (before and after
% %       the metric rectification) with the chosen lines printed on it.
% %       Compute also the angles between the pair of lines before and after
% %       rectification.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 4. OPTIONAL: Metric Rectification in a single step
% % Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)
% 
% %% 5. OPTIONAL: Affine Rectification of the left facade of image 0000
% 
% %% 6. OPTIONAL: Metric Rectification of the left facade of image 0000
% 
% %% 7. OPTIONAL: Affine Rectification of the left facade of image 0001
% 
% %% 8. OPTIONAL: Metric Rectification of the left facade of image 0001
% 

