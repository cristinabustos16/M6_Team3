% %% 2. Affine Rectification

% I = imread('Data/0000_s.png');
% A = load('Data/0000_s_info_lines.txt');
% affineRectification(I, A, 424, 240, 712, 565);

% The function "affineRectification" that gets as input an image I, the lines 
% information, and the indices of the lines to use, and computes and shows an
% affine srectification on the image.

function [transformedImage] = affineRectification(I, A, i1, i2, i3, i4)

    % indices of lines
    p1 = [A(i1,1) A(i1,2) 1]';
    p2 = [A(i1,3) A(i1,4) 1]';
    p3 = [A(i2,1) A(i2,2) 1]';
    p4 = [A(i2,3) A(i2,4) 1]';
    p5 = [A(i3,1) A(i3,2) 1]';
    p6 = [A(i3,3) A(i3,4) 1]';
    p7 = [A(i4,1) A(i4,2) 1]';
    p8 = [A(i4,3) A(i4,4) 1]';

    % ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
    l1 = cross(p1,p2);
    l2 = cross(p3,p4);
    l3 = cross(p5,p6);
    l4 = cross(p7,p8);

    % show the chosen lines in the image
    figure;imshow(I);
    hold on;
    t=1:0.1:1000;
    plot(t, -(l1(1)*t + l1(3)) / l1(2), 'r');
    plot(t, -(l2(1)*t + l2(3)) / l2(2), 'g');
    plot(t, -(l3(1)*t + l3(3)) / l3(2), 'b');
    plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

    % % ToDo: compute the homography that affinely rectifies the image

    % compute the vanishing points, lecture 2 slide 16
    v1 = cross(l1, l2);
    v2 = cross(l3, l4);

    % compute the line at infinity
    l = cross(v1, v2);
    % normalize line at infinity
    l = l / norm(l);

    % compute H based on the line at infinity, lecture2 slide 12
    H = [1 0 0; 0 1 0; l'];


    I2 = apply_H(permute(I,[2 1 3]), H);
    I2 = permute(I2,[2 1 3]);
    figure; imshow(uint8(I2));

    % ToDo: compute the transformed lines lr1, lr2, lr3, lr4

    % Transpose of the inverse of H:
    H_inv_tras = (inv(H))';

    lr1 = H_inv_tras * l1;
    lr2 = H_inv_tras * l2;
    lr3 = H_inv_tras * l3;
    lr4 = H_inv_tras * l4;

    % show the transformed lines in the transformed image
    figure;imshow(uint8(I2));
    hold on;
    t=1:0.1:1000;
    plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'r');
    plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'g');
    plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'b');
    plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

    % ToDo: to evaluate the results, compute the angle between the different pair 
    % of lines before and after the image transformation

    % Lines before transformation in cartesian coordinates:
    l1c = [l1(1) / l1(3), l1(2) / l1(3)]';
    l2c = [l2(1) / l2(3), l2(2) / l2(3)]';
    l3c = [l3(1) / l3(3), l3(2) / l3(3)]';
    l4c = [l4(1) / l4(3), l4(2) / l4(3)]';

    % Angles of lines before transformation:
    ang_l1_l2_before = acos((l1c'*l2c)/sqrt((l1c'*l1c)*(l2c'*l2c)));
    ang_l3_l4_before = acos((l3c'*l4c)/sqrt((l3c'*l3c)*(l4c'*l4c)));
    fprintf('Angle of l1 and l2 before transformation = %f ?\n', ang_l1_l2_before)
    fprintf('Angle of l3 and l4 before transformation = %f ?\n', ang_l3_l4_before)

    % Lines after transformation in cartesian coordinates:
    l1rc = [lr1(1) / lr1(3), lr1(2) / lr1(3)]';
    l2rc = [lr2(1) / lr2(3), lr2(2) / lr2(3)]';
    l3rc = [lr3(1) / lr3(3), lr3(2) / lr3(3)]';
    l4rc = [lr4(1) / lr4(3), lr4(2) / lr4(3)]';

    % Angles of lines after transformation:
    ang_l1_l2_after = acos((lr1c'*lr2c)/sqrt((lr1c'*lr1c)*(lr2c'*lr2c)));
    ang_l3_l4_after = acos((lr3c'*lr4c)/sqrt((lr3c'*lr3c)*(lr4c'*lr4c)));
    fprintf('Angle of l1 and l2 after transformation = %f ?\n', ang_l1_l2_after)
    fprintf('Angle of l3 and l4 after transformation = %f ?\n', ang_l3_l4_after)
end