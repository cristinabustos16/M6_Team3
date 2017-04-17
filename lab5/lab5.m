%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 5: Reconstruction from uncalibrated viewas


addpath('../lab2/sift'); % ToDo: change 'sift' to the correct path where you have the sift functions
addpath(genpath('./vanishing_points_v0.4/'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Create synthetic data

% 3D points
X = [];
Xi = [0 0 10]'; % middle facade
X = [X Xi];
Xi = [0 20 10]';
X = [X Xi];
Xi = [0 20 0]';
X = [X Xi];
Xi = [0 0 0]';
X = [X Xi];
Xi = [10 20 10]'; % right facade
X = [X Xi];
Xi = [10 20 0]';
X = [X Xi];
Xi = [30 0 10]';  % left facade
X = [X Xi];
Xi = [30 0 0]'; 
X = [X Xi];
Xi = [0 5 8]';   % middle squared window
X = [X Xi];
Xi = [0 8 8]';
X = [X Xi];
Xi = [0 5 5]';
X = [X Xi];
Xi = [0 8 5]';
X = [X Xi];
Xi = [0 12 5]';  % middle rectangular window
X = [X Xi];
Xi = [0 18 5]';
X = [X Xi];
Xi = [0 12 2]';
X = [X Xi];
Xi = [0 18 2]';
X = [X Xi];
Xi = [3 20 7]';   % left squared window
X = [X Xi];
Xi = [7 20 7]';
X = [X Xi];
Xi = [3 20 3]';
X = [X Xi];
Xi = [7 20 3]';
X = [X Xi];
Xi = [5 0 7]';   % right rectangular window
X = [X Xi];
Xi = [25 0 7]';
X = [X Xi];
Xi = [5 0 3]';
X = [X Xi];
Xi = [25 0 3]';
X = [X Xi];


% cameras

K = [709 0 450; 0 709 300; 0 0 1];
Rz = [cos(0.88*pi/2) -sin(0.88*pi/2) 0; sin(0.88*pi/2) cos(0.88*pi/2) 0; 0 0 1];
Ry = [cos(0.88*pi/2) 0 sin(0.88*pi/2); 0 1 0; -sin(0.88*pi/2) 0 cos(0.88*pi/2)];
R1 = Rz*Ry;
t1 = -R1*[40; 10; 5];

Rz = [cos(0.8*pi/2) -sin(0.8*pi/2) 0; sin(0.8*pi/2) cos(0.8*pi/2) 0; 0 0 1];
Ry = [cos(0.88*pi/2) 0 sin(0.88*pi/2); 0 1 0; -sin(0.88*pi/2) 0 cos(0.88*pi/2)];
Rx = [1 0 0; 0 cos(-0.15) -sin(-0.15); 0 sin(-0.15) cos(-0.15)];
R2 = Rx*Rz*Ry;
t2 = -R2*[45; 15; 5];

P1 = zeros(3,4);
P1(1:3,1:3) = R1;
P1(:,4) = t1;
P1 = K*P1;

P2 = zeros(3,4);
P2(1:3,1:3) = R2;
P2(:,4) = t2;
P2 = K*P2;

w = 900;
h = 600;

% visualize as point cloud
figure; hold on;
plot_camera2(P1,w,h);
plot_camera2(P2,w,h);
for i = 1:length(X)
    scatter3(X(1,i), X(2,i), X(3,i), 5^2, [0.5 0.5 0.5], 'filled');
end;
axis equal;
axis vis3d;

%% visualize as lines
figure;
hold on;
X1 = X(:,1); X2 = X(:,2); X3 = X(:,3); X4 = X(:,4);
plot3([X1(1) X2(1)], [X1(2) X2(2)], [X1(3) X2(3)]);
plot3([X3(1) X4(1)], [X3(2) X4(2)], [X3(3) X4(3)]);
X5 = X(:,5); X6 = X(:,6); X7 = X2; X8 = X3;
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,7); X6 = X(:,8); X7 = X1; X8 = X4;
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,9); X6 = X(:,10); X7 = X(:,11); X8 = X(:,12);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,13); X6 = X(:,14); X7 = X(:,15); X8 = X(:,16);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,17); X6 = X(:,18); X7 = X(:,19); X8 = X(:,20);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,21); X6 = X(:,22); X7 = X(:,23); X8 = X(:,24);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
axis vis3d
axis equal
plot_camera2(P1,w,h);
plot_camera2(P2,w,h);

%% Create homogeneous coordinates

% homogeneous 3D coordinates
Xh=[X; ones(1,length(X))];

% homogeneous 2D coordinates
x1 = P1*Xh;
x1(1,:) = x1(1,:)./x1(3,:);
x1(2,:) = x1(2,:)./x1(3,:);
x1(3,:) = x1(3,:)./x1(3,:);
x2 = P2*Xh;
x2(1,:) = x2(1,:)./x2(3,:);
x2(2,:) = x2(2,:)./x2(3,:);
x2(3,:) = x2(3,:)./x2(3,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Projective reconstruction (synthetic data)

% ToDo: create the function 'factorization_method' that computes a
% projective reconstruction with the factorization method of Sturm and
% Triggs '1996
% This function returns an estimate of:
%       Pproj: 3*Npoints x Ncam matrix containing the camera matrices
%       Xproj: 4 x Npoints matrix of homogeneous coordinates of 3D points
% 
% As a convergence criterion you may compute the Euclidean
% distance (d) between data points and projected points in both images 
% and stop when (abs(d - d_old)/d) < 0.1 where d_old is the distance
% in the previous iteration.

[Pproj, Xproj] = factorization_method({x1, x2});

%% Check projected points (estimated and data points)

for i=1:2
    x_proj{i} = euclid(Pproj(3*i-2:3*i, :)*Xproj);
end
x_d{1} = euclid(P1*Xh);
x_d{2} = euclid(P2*Xh);

% image 1
figure;
hold on
plot(x_d{1}(1,:),x_d{1}(2,:),'r*');
plot(x_proj{1}(1,:),x_proj{1}(2,:),'bo');
axis equal

% image 2
figure;
hold on
plot(x_d{2}(1,:),x_d{2}(2,:),'r*');
plot(x_proj{2}(1,:),x_proj{2}(2,:),'bo');


%% Visualize projective reconstruction
Xaux(1,:) = Xproj(1,:)./Xproj(4,:);
Xaux(2,:) = Xproj(2,:)./Xproj(4,:);
Xaux(3,:) = Xproj(3,:)./Xproj(4,:);
X=Xaux;

figure;
title('task 1: projective reconstruction');
hold on;
X1 = X(:,1); X2 = X(:,2); X3 = X(:,3); X4 = X(:,4);
plot3([X1(1) X2(1)], [X1(2) X2(2)], [X1(3) X2(3)]);
plot3([X3(1) X4(1)], [X3(2) X4(2)], [X3(3) X4(3)]);
X5 = X(:,5); X6 = X(:,6); X7 = X2; X8 = X3;
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,7); X6 = X(:,8); X7 = X1; X8 = X4;
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,9); X6 = X(:,10); X7 = X(:,11); X8 = X(:,12);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,13); X6 = X(:,14); X7 = X(:,15); X8 = X(:,16);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,17); X6 = X(:,18); X7 = X(:,19); X8 = X(:,20);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,21); X6 = X(:,22); X7 = X(:,23); X8 = X(:,24);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
axis vis3d
axis equal



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine reconstruction (synthetic data)

% ToDo: create the function 'vanishing_point' that computes the vanishing
% point formed by the line that joins points xo1 and xf1 and the line 
% that joins points x02 and xf2
% 
% lpmayos note: lecture 9a
    % There are different methods for autocalibration (process of determining 
    % the internal camera parameters directly from multiple uncalibrated images). 
    % We will see the stratified methods which involve two steps:
    %   ??? Affine reconstruction: it finds the plane at infinity
    %   ??? Metric reconstruction: it finds the image of the absolute conic ?? = K ???T K ???1 (and thus the internal parameters)
    % 
    % We need to find the homography:
    %       Ha???p = [I 0; pT 1]
    %       The essence of the affine reconstruction is then to locate the plane at 
    %       infinity in the projective reconstruction frame.
    %       By applying Ha???p we will map p to ????? = (0,0,0,1)T and we will recover te paralellism.
    % How do we identify p?  (slide 22)
    %       As every plane, the plane at infinity is determined by three points on it.
    % since ??TX_i = 0, i = 1,2,3, we have: 
    %       X_1^T
    %       X_2^T  ?? = 0 --> A?? = 0  where A is a 3x4 matrix
    %       X_3^T
    % If the three points are in general position (not on the same line), they
    % provide linearly indep. eq. and the matrix they form is rank 3.
    % ?????? Then ?? is obtained uniquely (up to scale) as the 1-dimensional right null space of A.

%
% [v1] = vanishing_point(xo1, xf1, xo2, xf2)

% Compute the vanishing points in each image
v1 = vanishing_point(x1(:,21),x1(:,22),x1(:,23),x1(:,24));
v2 = vanishing_point(x1(:,21),x1(:,23),x1(:,22),x1(:,24));
v3 = vanishing_point(x1(:,1),x1(:,2),x1(:,4),x1(:,3));

v3(3) = -1.8190e-11; 

v1p = vanishing_point(x2(:,21),x2(:,22),x2(:,23),x2(:,24));
v2p = vanishing_point(x2(:,21),x2(:,23),x2(:,22),x2(:,24));
v3p = vanishing_point(x2(:,1),x2(:,2),x2(:,4),x2(:,3));

% ToDo: use the vanishing points to compute the matrix Hp that 
%       upgrades the projective reconstruction to an affine reconstruction

% .............................
p1 = Pproj(1:3,:);  % camera 1
p2 = Pproj(4:6,:);  % camera 2

A = [triangulate(euclid(v1), euclid(v1p), p1, p2, [w h])'; ...
     triangulate(euclid(v2), euclid(v2p), p1, p2, [w h])'; ...
     triangulate(euclid(v3), euclid(v3p), p1, p2, [w h])'];

% lpmayos note
% Then ?? is obtained uniquely (up to scale) as the 1-dimensional right null space of A.
% We can compute the right null space of A with the svd (ref: https://cseweb.ucsd.edu/classes/wi15/cse252B-a/nullspace.pdf)
% ... i-th column of V are the corresponding right singular vector of A.
% ... The rank r of A is the number of nonzero singular values. 
% ... The (right) null space of A is the columns of V corresponding to singular values equal to zero.

rightNullSpace = null(A);  % Z = null(A) is an orthonormal basis for the null space of A obtained from the singular value decomposition.
H = [euclid(rightNullSpace)' 1];  % we could just divide each element by H(4)
% equivalent to [U, D, V] = svd(A); H = V(:, end);
Hp = eye(4);
Hp(end,:) = H;

% .............................

%% check results

Xa = euclid(Hp*Xproj);
figure;
title('task 2 affine reconstruction');
hold on;
X1 = Xa(:,1); X2 = Xa(:,2); X3 = Xa(:,3); X4 = Xa(:,4);
plot3([X1(1) X2(1)], [X1(2) X2(2)], [X1(3) X2(3)]);
plot3([X3(1) X4(1)], [X3(2) X4(2)], [X3(3) X4(3)]);
X5 = Xa(:,5); X6 = Xa(:,6); X7 = X2; X8 = X3;
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,7); X6 = Xa(:,8); X7 = X1; X8 = X4;
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,9); X6 = Xa(:,10); X7 = Xa(:,11); X8 = Xa(:,12);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,13); X6 = Xa(:,14); X7 = Xa(:,15); X8 = Xa(:,16);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,17); X6 = Xa(:,18); X7 = Xa(:,19); X8 = Xa(:,20);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,21); X6 = Xa(:,22); X7 = Xa(:,23); X8 = Xa(:,24);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
axis vis3d
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric reconstruction (synthetic data)

% ToDo: compute the matrix Ha that 
%       upgrades the projective reconstruction to an affine reconstruction
% Use the following vanishing points given by three pair of orthogonal lines
% and assume that the skew factor is zero and that pixels are square

v1 = vanishing_point(x1(:,2),x1(:,5),x1(:,3),x1(:,6));
v2 = vanishing_point(x1(:,1),x1(:,2),x1(:,3),x1(:,4));
v3 = vanishing_point(x1(:,1),x1(:,4),x1(:,2),x1(:,3));


% lpmayos note: see nice intro to metric reconstruction on book page 272 (pdf 290)

% book pdf page 295 algorithm 10.1



% lpmayos note: lecture 9a, slide 29 (book page 478 / pdf page 496)
% lpmayos note: see algorithm 19.2 on book page 479 / pdf page 497
%       Metric rectification: Determine the camera matrix K from the Cholesky 
%       decomposition ?? = (KK^T)^???1 . Then a metric reconstruction is obtained
%       as {Pi H_P H_A, (H_P H_A )???1 X_j } with H_A = [K 0; 0^T 1]


% note: notation xxa, a added to differentiate matrix elements from vanishing points
u1a = v1(1); u2a = v1(2); u3a = v1(3);
v1a = v2(1); v2a = v2(2); v3a = v2(3);
z1a = v3(1); z2a = v3(2); z3a = v3(3);

A = [ u1a*v1a u1a*v2a+u2a*v1a u1a*v3a+u3a*v1a u2a*v2a u2a*v3a+u3a*v2a u3a*v3a; ...
      u1a*z1a u1a*z2a+u2a*z1a u1a*z3a+u3a*z1a u2a*z2a u2a*z3a+u3a*z2a u3a*z3a; ...
      v1a*z1a v1a*z2a+v2a*z1a v1a*z3a+v3a*z1a v2a*z2a v2a*z3a+v3a*z2a v3a*z3a; ...
      0       1               0               0       0               0; ...
      1       0               0               -1      0               0 ];

% book 496: Algorithm 19.2 comes down to solving a homogeneous set of equations
% of the form Ac = 0, where c represents w arranged as a 6-vector.
% In matrix form: A??V = 0, where ??V = (??11,??12,??13,??22,??23,??33)T

% lpmayos: The solution ??_V is the null vector of A (slides 9a slide 35).
wv = null(A);
wv = wv(:,end);

% lpmayos: slides 9a slide 34
omega = [ wv(1) wv(2) wv(3); ...
      wv(2) wv(4) wv(5); ...
      wv(3) wv(5) wv(6) ];

% forcing w to be positive definite  
% [vec,val]=eig(w);
% val(val<0)=eps;
% w=vec*val*vec';

% the camera is known to have zero skew, then w12 = 0
% pixels are square, that is, zero skew and alpha_x = alpha_y , then: w11 = w22
% w = [ wv(2,2) 0       wv(1,3); ...
%       0       wv(2,2) wv(2,3); ...
%       wv(1,3) wv(2,3) wv(3,3) ];

% book alpg. 10.1 page pdf 295
% camera in the affine reconstruction for which w is computed
p1_a = p1 * inv(Hp);
% p1_a = p1 * Hp;
% M is the first 3 ?? 3 submatrix
M = p1_a(1:3, 1:3);
% A is obtained by Cholesky factorization from the equation AA^T = (M^T w M)^???1


A = chol(inv(M' * omega * M));
% Homography to upgrade the affine reconstruction to a metric reconstruction
Ha = [inv(A), [0, 0, 0]'; 0, 0, 0, 1];


%% check results

Xa = euclid(Ha*Hp*Xproj);
figure;
hold on;
X1 = Xa(:,1); X2 = Xa(:,2); X3 = Xa(:,3); X4 = Xa(:,4);
plot3([X1(1) X2(1)], [X1(2) X2(2)], [X1(3) X2(3)]);
plot3([X3(1) X4(1)], [X3(2) X4(2)], [X3(3) X4(3)]);
X5 = Xa(:,5); X6 = Xa(:,6); X7 = X2; X8 = X3;
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,7); X6 = Xa(:,8); X7 = X1; X8 = X4;
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,9); X6 = Xa(:,10); X7 = Xa(:,11); X8 = Xa(:,12);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,13); X6 = Xa(:,14); X7 = Xa(:,15); X8 = Xa(:,16);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,17); X6 = Xa(:,18); X7 = Xa(:,19); X8 = Xa(:,20);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,21); X6 = Xa(:,22); X7 = Xa(:,23); X8 = Xa(:,24);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
axis vis3d
axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Projective reconstruction (real data)

%% read images
Irgb{1} = double(imread('Data/0000_s.png'))/255;
Irgb{2} = double(imread('Data/0001_s.png'))/255;

I{1} = sum(Irgb{1}, 3) / 3; 
I{2} = sum(Irgb{2}, 3) / 3;

Ncam = length(I);

% ToDo: compute a projective reconstruction using the factorization method

% ToDo: show the data points (image correspondences) and the projected
% points (of the reconstructed 3D points) in images 1 and 2. Reuse the code
% in section 'Check projected points' (synthetic experiment).

% Compute sift keypoints and matches.
points = cell(2,1);
descriptors = cell(2,1);
for i = 1:Ncam
    [points{i}, descriptors{i}] = sift(I{i}, 'Threshold', 0.01);
    points{i} = points{i}(1:2,:);
end

matches = siftmatch(descriptors{1}, descriptors{2});
% Plot matches.
figure();
plotmatches(I{1}, I{2}, points{1}, points{2}, matches, 'Stacking', 'v');

x1 = points{1}(:,matches(1,:));
x2 = points{2}(:,matches(2,:));
[F, inliers] = ransac_fundamental_matrix(homog(x1), homog(x2), 2);

% Plot inliers.
inlier_matches = matches(:, inliers);
figure;
plotmatches(I{1}, I{2}, points{1}, points{2}, inlier_matches, 'Stacking', 'v');

x1 = points{1}(:, inlier_matches(1, :));
x2 = points{2}(:, inlier_matches(2, :));

x1 = homog(x1);
x2 = homog(x2);

[Pproj, Xproj] = factorization_method({x1, x2});
% Check projected points (estimated and data points)
for i=1:2
    x_proj{i} = euclid(Pproj(3*i-2:3*i, :)*Xproj);
end
x_d{1} = euclid(P1*Xh);
x_d{2} = euclid(P2*Xh);

% image 1
figure;
hold on
plot(x_d{1}(1,:),x_d{1}(2,:),'r*');
plot(x_proj{1}(1,:),x_proj{1}(2,:),'bo');
axis equal

% image 2
figure;
hold on
plot(x_d{2}(1,:),x_d{2}(2,:),'r*');
plot(x_proj{2}(1,:),x_proj{2}(2,:),'bo');

%3D visualisation
% X(1,:) = Xproj(1,:)./Xproj(4,:);
% X(2,:) = Xproj(2,:)./Xproj(4,:);
% X(3,:) = Xproj(3,:)./Xproj(4,:);
% figure; hold on;
% for i = 1:length(X)
%     scatter3(X(1,i), X(2,i), X(3,i), 15, [0.2 0.9 0.4], 'filled');
% end;
% axis([-200 200 -200 200 -200 200]);
% axis vis3d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Affine reconstruction (real data)

% ToDo: compute the matrix Hp that updates the projective reconstruction
% to an affine one
%
% You may use the vanishing points given by function 'detect_vps' that 
% implements the method presented in Lezama et al. CVPR 2014
% (http://dev.ipol.im/~jlezama/vanishing_points/)

% This is an example on how to obtain the vanishing points (VPs) from three
% orthogonal lines in image 1

img_in1 =  'Data/0000_s.png'; % input image
folder_out = '.'; % output folder
manhattan = 1;
acceleration = 0;
focal_ratio = 1;
params.PRINT = 1;
params.PLOT = 1;
[horizon1, VPs1] = detect_vps(img_in, folder_out, manhattan, acceleration, focal_ratio, params);

img_in2 =  'Data/0001_s.png'; % input image
[horizon2, VPs2] = detect_vps(img_in2, folder_out, manhattan, acceleration, focal_ratio, params);

p1 = Pproj(1:3,:);  % camera 1
p2 = Pproj(4:6,:);  % camera 2

[w,h] = size(I{1});

A = [triangulate(VPs1(:,1), VPs2(:,1), p1, p2, [w h])'; ...
     triangulate(VPs1(:,2), VPs2(:,2), p1, p2, [w h])'; ...
     triangulate(VPs1(:,3), VPs2(:,3), p1, p2, [w h])'];

rightNullSpace = null(A);  % Z = null(A) is an orthonormal basis for the null space of A obtained from the singular value decomposition.
H = [euclid(rightNullSpace)' 1];  % we could just divide each element by H(4)
% equivalent to [U, D, V] = svd(A); H = V(:, end);
Hp = eye(4);
Hp(end,:) = H;


%% Visualize the result

% x1m are the data points in image 1
% Xm are the reconstructed 3D points (projective reconstruction)
x1m = euclid(x1);
Xm = Xproj;

r = interp2(double(Irgb{1}(:,:,1)), x1m(1,:), x1m(2,:));
g = interp2(double(Irgb{1}(:,:,2)), x1m(1,:), x1m(2,:));
b = interp2(double(Irgb{1}(:,:,3)), x1m(1,:), x1m(2,:));
Xe = euclid(Hp*Xm);
figure; hold on;
[w,h] = size(I{1});
for i = 1:length(Xe)
    scatter3(Xe(1,i), Xe(2,i), Xe(3,i), 2^2, [r(i) g(i) b(i)], 'filled');
end;
axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Metric reconstruction (real data)

% ToDo: compute the matrix Ha that updates the affine reconstruction
% to a metric one and visualize the result in 3D as in the previous section
v1 = homog(VPs1(:,1));
v2 = homog(VPs1(:,2));
v3 = homog(VPs1(:,3));

% note: notation xxa, a added to differentiate matrix elements from vanishing points
u1a = v1(1); u2a = v1(2); u3a = v1(3);
v1a = v2(1); v2a = v2(2); v3a = v2(3);
z1a = v3(1); z2a = v3(2); z3a = v3(3);

A = [ u1a*v1a u1a*v2a+u2a*v1a u1a*v3a+u3a*v1a u2a*v2a u2a*v3a+u3a*v2a u3a*v3a; ...
      u1a*z1a u1a*z2a+u2a*z1a u1a*z3a+u3a*z1a u2a*z2a u2a*z3a+u3a*z2a u3a*z3a; ...
      v1a*z1a v1a*z2a+v2a*z1a v1a*z3a+v3a*z1a v2a*z2a v2a*z3a+v3a*z2a v3a*z3a; ...
      0       1               0               0       0               0; ...
      1       0               0               -1      0               0 ];
  
wv = null(A);
wv = wv(:,end);

omega = [ wv(1) wv(2) wv(3); ...
          wv(2) wv(4) wv(5); ...
          wv(3) wv(5) wv(6) ];

% forcing w to be positive definite  
% [vec,val]=eig(w);
% val(val<0)=eps;
% w=vec*val*vec';

% book alpg. 10.1 page pdf 295
% camera in the affine reconstruction for which w is computed
p1_a = p1 * inv(Hp);

% p1_a = p1 * Hp;
% M is the first 3 ?? 3 submatrix

M = p1_a(1:3, 1:3);

% A is obtained by Cholesky factorization from the equation AA^T = (M^T w M)^???1
A = chol(inv(M' * omega * M));

% Homography to upgrade the affine reconstruction to a metric reconstruction
Ha = [inv(A), [0, 0, 0]'; 0, 0, 0, 1];

%% Visualize the result

% x1m are the data points in image 1
% Xm are the reconstructed 3D points (projective reconstruction)

r = interp2(double(Irgb{1}(:,:,1)), x1m(1,:), x1m(2,:));
g = interp2(double(Irgb{1}(:,:,2)), x1m(1,:), x1m(2,:));
b = interp2(double(Irgb{1}(:,:,3)), x1m(1,:), x1m(2,:));
Xe = euclid(Ha*Hp*Xm);
figure; hold on;
[w,h] = size(I{1});
for i = 1:length(Xe)
    scatter3(Xe(1,i), Xe(2,i), Xe(3,i), 2^2, [r(i) g(i) b(i)], 'filled');
end;
axis equal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. OPTIONAL: Projective reconstruction from two views

% ToDo: compute a projective reconstruction from the same two views 
% by computing two possible projection matrices from the fundamental matrix
% and one of the epipoles.
% Then update the reconstruction to affine and metric as before (reuse the code).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8. OPTIONAL: Projective reconstruction from more than two views

% ToDo: extend the function that computes the projective reconstruction 
% with the factorization method to the case of three views. You may use 
% the additional image '0002_s.png'
% Then update the reconstruction to affine and metric.
%
% Any other improvement you may icorporate (add a 4th view,
% incorporate new 3D points by triangulation, incorporate new views
% by resectioning, better visualization of the result with another
% software, ...)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 9. OPTIONAL: Any other improvement you may icorporate 
%
%  (add a 4th view, incorporate new 3D points by triangulation, 
%   apply any kind of processing on the point cloud, ...)
