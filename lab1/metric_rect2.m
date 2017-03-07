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
i = 534;
p9 = [A(i,1) A(i,2) 1]';
p10 = [A(i,3) A(i,4) 1]';
i = 576;
p11 = [A(i,1) A(i,2) 1]';
p12 = [A(i,3) A(i,4) 1]';
i = 367;
p13 = [A(i,1) A(i,2) 1]';
p14 = [A(i,3) A(i,4) 1]';
i = 227;
p15 = [A(i,1) A(i,2) 1]';
p16 = [A(i,3) A(i,4) 1]';

% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = cross(p1,p2);
l2 = cross(p3,p4);
l3 = cross(p5,p6);
l4 = cross(p7,p8);
l5 = cross(p9,p10);
l6 = cross(p11,p12);
l7 = cross(p13,p14);
l8 = cross(p15,p16);

B = zeros(5,6);
l = l1;
m = l3;
B(1,:) = [l(1)*m(1), (l(1)*m(2) + l(2)*m(2))/2, l(2)*m(2), (l(1)*m(3) + l(3)*m(2))/2, (l(2)*m(3) + l(3)*m(2))/2, l(3)*m(3)];
l = l2;
m = l4;
B(2,:) = [l(1)*m(1), (l(1)*m(2) + l(2)*m(2))/2, l(2)*m(2), (l(1)*m(3) + l(3)*m(2))/2, (l(2)*m(3) + l(3)*m(2))/2, l(3)*m(3)];
l = l5;
m = l7;
B(3,:) = [l(1)*m(1), (l(1)*m(2) + l(2)*m(2))/2, l(2)*m(2), (l(1)*m(3) + l(3)*m(2))/2, (l(2)*m(3) + l(3)*m(2))/2, l(3)*m(3)];
l = l6;
m = l8;
B(4,:) = [l(1)*m(1), (l(1)*m(2) + l(2)*m(2))/2, l(2)*m(2), (l(1)*m(3) + l(3)*m(2))/2, (l(2)*m(3) + l(3)*m(2))/2, l(3)*m(3)];
l = l1;
m = l5;
B(5,:) = [l(1)*m(1), (l(1)*m(2) + l(2)*m(2))/2, l(2)*m(2), (l(1)*m(3) + l(3)*m(2))/2, (l(2)*m(3) + l(3)*m(2))/2, l(3)*m(3)];

[U, S, V] = svd(B, 'econ');
c = V(:,end); % I hope this is the null vector of B.
C_inf_star = [c(1)  , c(2)/2, c(4)/2;
              c(2)/2, c(3)  , c(5)/2;
              c(4)/2, c(5)/2, c(6)];
          
[U, S, V] = svd(C_inf_star);

H = U;

% Transform the lines:
l1r = inv(H)' * l1;
l2r = inv(H)' * l2;
l3r = inv(H)' * l3;

compute_angle(l1r, l2r) % this should be 0
compute_angle(l1r, l3r) % this should be 90



