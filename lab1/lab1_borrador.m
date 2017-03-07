%%%%%%%%%%%%%%%%%%
% choose the image points
I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% indices of lines
i = 424;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 240;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 712;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 565;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = cross(p1,p2);
l2 = cross(p3,p4);
l3 = cross(p5,p6);
l4 = cross(p7,p8);

% show the chosen lines in the image
% figure;imshow(I);
% hold on;
% t=1:0.1:1000;
% plot(t, -(l1(1)*t + l1(3)) / l1(2), 'r');
% plot(t, -(l2(1)*t + l2(3)) / l2(2), 'g');
% plot(t, -(l3(1)*t + l3(3)) / l3(2), 'b');
% plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% % ToDo: compute the homography that affinely rectifies the image

% compute the vanishing points, lecture 2 slide 16
v1 = cross(l1, l2);
v2 = cross(l3, l4);

% compute the line at infinity
linf = cross(v1, v2);
% normalize line at infinity
linf = linf / max (linf);

% compute H based on the line at infinity, lecture2 slide 12
H_p = [1 0 0; 0 1 0; linf'];


% ToDo: compute the transformed lines lr1, lr2, lr3, lr4

% Transpose of the inverse of H:
H_p_inv_tras = (inv(H_p))';

l1_a = H_p_inv_tras * l1;
l2_a = H_p_inv_tras * l2;
l3_a = H_p_inv_tras * l3;
l4_a = H_p_inv_tras * l4;

% Angles of lines before transformation:
ang_l1_l2_before = compute_angle(l1, l2);
ang_l3_l4_before = compute_angle(l3, l4);
fprintf('Angle of l1 and l2 before transformation = %f �\n', ang_l1_l2_before)
fprintf('Angle of l3 and l4 before transformation = %f �\n', ang_l3_l4_before)

% Angles of lines after transformation:
ang_l1_l2_after = compute_angle(l1_a, l2_a);
ang_l3_l4_after = compute_angle(l3_a, l4_a);
fprintf('Angle of l1 and l2 after transformation = %f �\n', ang_l1_l2_after)
fprintf('Angle of l3 and l4 after transformation = %f �\n', ang_l3_l4_after)


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
s = V(:,2); % I hope this is the null vector of B.
S = [s(1), s(2);
     s(2), s(3)];

% Compute the upper triangular matrix using the Cholesky factorization:
A = chol(S, 'lower');
% A = A / sqrt(det(A));

% Matrix for the metric rectification:
H_a = [A(1,1), A(1,2), 0;
       A(2,1), A(2,2), 0;
       0     , 0     , 1];

% H = H_p * H_a;

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
% l1_s = inv(H)' * l1;
% l2_s = inv(H)' * l2;
% l3_s = inv(H)' * l3;
% l4_s = inv(H)' * l4;

l1_s = inv(H_a)' * l1_a;
l2_s = inv(H_a)' * l2_a;
l3_s = inv(H_a)' * l3_a;
l4_s = inv(H_a)' * l4_a;

compute_angle(l1_s, l3_s)
compute_angle(l2_s, l4_s)

% Angles of lines before transformation:
ang_l1_l3_before = compute_angle(l1_a, l3_a);
ang_l2_l4_before = compute_angle(l2_a, l4_a);
fprintf('Angle of l1 and l3 before metric rectification = %f �\n', ang_l1_l3_before)
fprintf('Angle of l2 and l4 before metric rectification = %f �\n', ang_l2_l4_before)

% Angles of lines after transformation:
ang_l1_l3_after = compute_angle(l1_s, l3_s);
ang_l2_l4_after = compute_angle(l2_s, l4_s);
fprintf('Angle of l1 and l3 after metric rectification = %f �\n', ang_l1_l3_after)
fprintf('Angle of l2 and l4 after metric rectification = %f �\n', ang_l2_l4_after)




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

