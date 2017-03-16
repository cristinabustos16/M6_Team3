%%%% DLT
% Aplicamos el algoritmo DLT para rectificar la imagen de la semana 1.
% Usamos 8 puntos, de modo que tenemos un sistema sobredeterminado.

clearvars;
close all

addpath('./../lab1/')

% Carga datos de la semana 1:
I=imread('./../lab1/Data/0000_s.png');
A = load('./../lab1/Data/0000_s_info_lines.txt');

% Indices de las lineas:
i1 = 424;
i2 = 240;
i3 = 712;
i4 = 565;
i5 = 227;
i6 = 367;
i7 = 534;
i8 = 576;

% Linea 1:
p1 = [A(i1,1) A(i1,2) 1]';
p2 = [A(i1,3) A(i1,4) 1]';
l1 = cross(p1,p2);

% Linea 2:
p3 = [A(i2,1) A(i2,2) 1]';
p4 = [A(i2,3) A(i2,4) 1]';
l2 = cross(p3,p4);

% Linea 3:
p5 = [A(i3,1) A(i3,2) 1]';
p6 = [A(i3,3) A(i3,4) 1]';
l3 = cross(p5,p6);

% Linea 4:
p7 = [A(i4,1) A(i4,2) 1]';
p8 = [A(i4,3) A(i4,4) 1]';
l4 = cross(p7,p8);

% Linea 5:
p5 = [A(i5,1) A(i5,2) 1]';
p6 = [A(i5,3) A(i5,4) 1]';
l5 = cross(p5,p6);

% Linea 6:
p7 = [A(i6,1) A(i6,2) 1]';
p8 = [A(i6,3) A(i6,4) 1]';
l6 = cross(p7,p8);

% Linea 7:
p9 = [A(i7,1) A(i7,2) 1]';
p10 = [A(i7,3) A(i7,4) 1]';
l7 = cross(p9,p10);

% Linea 8:
p11 = [A(i8,1) A(i8,2) 1]';
p12 = [A(i8,3) A(i8,4) 1]';
l8 = cross(p11,p12);

% Esquinas de las ventanas (intersecciones entre las lineas):
% Ventana 1:
x1 = cross(l1, l3);
x2 = cross(l1, l4);
x3 = cross(l2, l3);
x4 = cross(l2, l4);
% Ventana 1:
x5 = cross(l5, l7);
x6 = cross(l5, l8);
x7 = cross(l6, l7);
x8 = cross(l6, l8);

% Destino:
% Ventana 1:
x1p = [100;100;1];
x2p = [150;100;1];
x3p = [100;150;1];
x4p = [150;150;1];
% Ventana 2:
x5p = [480;160;1];
x6p = [530;160;1];
x7p = [480;210;1];
x8p = [530;210;1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    DLT    %%%%

% Source points (from input image):
xsrc = [x1, x2, x3, x4, x5, x6, x7, x8];
% Destination points:
xdes = [x1p, x2p, x3p, x4p, x5p, x6p, x7p, x8p];

% Compute the homography relating the two sets of points:
H = homography2d(xsrc, xdes);

% Transformed lines:
l1p = inv(H)' * l1;
l2p = inv(H)' * l2;
l3p = inv(H)' * l3;
l4p = inv(H)' * l4;
l5p = inv(H)' * l5;
l6p = inv(H)' * l6;
l7p = inv(H)' * l7;
l8p = inv(H)' * l8;

% Check angles:
compute_angle(l1p, l2p) % 0
compute_angle(l3p, l4p) % 0
compute_angle(l1p, l3p) % 90
compute_angle(l2p, l4p) % 90
compute_angle(l5p, l6p) % 0
compute_angle(l7p, l8p) % 0
compute_angle(l5p, l7p) % 90
compute_angle(l6p, l8p) % 90

% Coordinates for the new image:
xmin = -500;
xmax = 1100;
ymin = -500;
ymax = 600;

% xmin = 0;
% xmax = size(I,2);
% ymin = 0;
% ymax = size(I,1);

% Corners of the new figure.
% Select them so the image fits in the figure.
corners = [xmin, xmax, ymin, ymax];

% Transformation to apply to points and lines in order to plot them in the
% transformed image:
Hfit = compute_Hfit(corners);

% Transform lines to fit in the new image:
l1fit = inv(Hfit)' * l1p;
l2fit = inv(Hfit)' * l2p;
l3fit = inv(Hfit)' * l3p;
l4fit = inv(Hfit)' * l4p;
l5fit = inv(Hfit)' * l5p;
l6fit = inv(Hfit)' * l6p;
l7fit = inv(Hfit)' * l7p;
l8fit = inv(Hfit)' * l8p;

% Transform points to fit in the new image:
x1fit = Hfit * (H * x1);
x2fit = Hfit * (H * x2);
x3fit = Hfit * (H * x3);
x4fit = Hfit * (H * x4);
x5fit = Hfit * (H * x5);
x6fit = Hfit * (H * x6);
x7fit = Hfit * (H * x7);
x8fit = Hfit * (H * x8);

% Transform image:
I2 = apply_H_v2(I, H, corners);

% Original image:
figure()
imshow(I)
hold on

% Window corners:
plot(x1(1)/x1(3), x1(2)/x1(3), '*r', 'LineWidth', 1)
plot(x2(1)/x2(3), x2(2)/x2(3), '*g', 'LineWidth', 1)
plot(x3(1)/x3(3), x3(2)/x3(3), '*b', 'LineWidth', 1)
plot(x4(1)/x4(3), x4(2)/x4(3), '*y', 'LineWidth', 1)
plot(x5(1)/x5(3), x5(2)/x5(3), '*r', 'LineWidth', 1)
plot(x6(1)/x6(3), x6(2)/x6(3), '*g', 'LineWidth', 1)
plot(x7(1)/x7(3), x7(2)/x7(3), '*b', 'LineWidth', 1)
plot(x8(1)/x8(3), x8(2)/x8(3), '*y', 'LineWidth', 1)

% Window lines:
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'r');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'g');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'b');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');
plot(t, -(l5(1)*t + l5(3)) / l5(2), 'r');
plot(t, -(l6(1)*t + l6(3)) / l6(2), 'g');
plot(t, -(l7(1)*t + l7(3)) / l7(2), 'b');
plot(t, -(l8(1)*t + l8(3)) / l8(2), 'y');

% Transformed image:
figure()
imshow(I2)
hold on

% Window corners:
plot(x1fit(1)/x1fit(3), x1fit(2)/x1fit(3), '*r', 'LineWidth', 1)
plot(x2fit(1)/x2fit(3), x2fit(2)/x2fit(3), '*g', 'LineWidth', 1)
plot(x3fit(1)/x3fit(3), x3fit(2)/x3fit(3), '*b', 'LineWidth', 1)
plot(x4fit(1)/x4fit(3), x4fit(2)/x4fit(3), '*y', 'LineWidth', 1)
plot(x5fit(1)/x5fit(3), x5fit(2)/x5fit(3), '*r', 'LineWidth', 1)
plot(x6fit(1)/x6fit(3), x6fit(2)/x6fit(3), '*g', 'LineWidth', 1)
plot(x7fit(1)/x7fit(3), x7fit(2)/x7fit(3), '*b', 'LineWidth', 1)
plot(x8fit(1)/x8fit(3), x8fit(2)/x8fit(3), '*y', 'LineWidth', 1)

% Window lines:
t=1:0.1:2000;
plot(t, -(l1fit(1)*t + l1fit(3)) / l1fit(2), 'r');
plot(t, -(l2fit(1)*t + l2fit(3)) / l2fit(2), 'g');
plot(t, -(l3fit(1)*t + l3fit(3)) / l3fit(2), 'b');
plot(t, -(l4fit(1)*t + l4fit(3)) / l4fit(2), 'y');
plot(t, -(l5fit(1)*t + l5fit(3)) / l5fit(2), 'r');
plot(t, -(l6fit(1)*t + l6fit(3)) / l6fit(2), 'g');
plot(t, -(l7fit(1)*t + l7fit(3)) / l7fit(2), 'b');
plot(t, -(l8fit(1)*t + l8fit(3)) / l8fit(2), 'y');

