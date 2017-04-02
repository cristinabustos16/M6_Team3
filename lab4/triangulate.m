% the function triangulate.m that performs a triangulation
%       with the homogeneous algebraic method (DLT)
%
%       The entries are (x1, x2, P1, P2, imsize), where:
%           - x1, and x2 are the Euclidean coordinates of two matching 
%             points in two different images.
%           - P1 and P2 are the two camera matrices
%           - imsize is a two-dimensional vector with the image size


function [ triang ] = triangulate( x1, x2, P1, P2, imsize )

% Normalization: Scaling and translation so that both pixel coordinates
% are in the interval [-1; 1].
nx = imsize(1);
ny = imsize(2);
H = [2/nx 0    -1; 
     0    2/ny -1; 
     0    0     1];

% We build A from x1 = Hx1, x2 = Hx2, P1 = HP1, and P2 = HP2.
x1 = H*homog(x1);
x2 = H*homog(x2);

P1 = H*P1;
P2 = H*P2;

x1_x = x1(1);
x1_y = x1(2);
x2_x = x2(1);
x2_y = x2(2);

A = [x1_x*P1(3,:) - P1(1,:); 
     x1_y*P1(3,:) - P1(2,:); 
     x2_x*P2(3,:) - P2(1,:); 
     x2_y*P2(3,:) - P2(2,:)];
 
[~,~,V] = svd(A);
triang = V(:,end);

end