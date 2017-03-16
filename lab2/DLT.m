%%%% DLT

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

% Esquinas de la ventana (intersecciones entre las lineas):
esq_sup_izq = cross(l1, l3);
esq_sup_der = cross(l1, l4);
esq_inf_izq = cross(l2, l3);
esq_inf_der = cross(l2, l4);

% Normalizamos los puntos:
esq_sup_izq = esq_sup_izq / esq_sup_izq(3);
esq_sup_der = esq_sup_der / esq_sup_der(3);
esq_inf_izq = esq_inf_izq / esq_inf_izq(3);
esq_inf_der = esq_inf_der / esq_inf_der(3);

%%%%%%%%%%%
%%% asignamos los puntos:

% Origen:
x1 = esq_sup_izq;
x2 = esq_sup_der;
x3 = esq_inf_izq;
x4 = esq_inf_der;

% Destino:
x1p = [100;100;1];
x2p = [150;100;1];
x3p = [100;150;1];
x4p = [150;150;1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    DLT    %%%%

% Matrices con las cuatro correspondencias:
A1 = [0, 0, 0     , -x1p(3) * x1', x1p(2) * x1';
      x1p(3) * x1', 0, 0, 0      , -x1p(1) * x1'];

A2 = [0, 0, 0     , -x2p(3) * x2', x2p(2) * x2';
      x2p(3) * x2', 0, 0, 0      , -x2p(1) * x2'];

A3 = [0, 0, 0     , -x3p(3) * x3', x3p(2) * x3';
      x3p(3) * x3', 0, 0, 0      , -x3p(1) * x3'];

A4 = [0, 0, 0     , -x4p(3) * x4', x4p(2) * x4';
      x4p(3) * x4', 0, 0, 0      , -x4p(1) * x4'];

% Ensamblamos las matrices:
A = [A1; A2; A3; A4];

% Solución del sistema, como vector nulo de A:
h = null(A);

% Construimos la transformación proyectiva:
H = [h(1), h(2), h(3);
     h(4), h(5), h(6);
     h(7), h(8), h(9)];

% Transformed lines:
l1p = inv(H)' * l1;
l2p = inv(H)' * l2;
l3p = inv(H)' * l3;
l4p = inv(H)' * l4;

% Check angles:
compute_angle(l1p, l2p) % 0
compute_angle(l3p, l4p) % 0
compute_angle(l1p, l3p) % 90
compute_angle(l2p, l4p) % 90

% Coordinates for the new image:
xmin = -500;
xmax = 1100;
ymin = -500;
ymax = 600;

% xmin = 0;
% xmax = size(I,2);
% ymin = 0;
% ymax = size(I,1);

corners = [xmin, xmax, ymin, ymax];

Hfit = compute_Hfit(corners);

% Transform lines to fit in the new image:
l1fit = inv(Hfit)' * l1p;
l2fit = inv(Hfit)' * l2p;
l3fit = inv(Hfit)' * l3p;
l4fit = inv(Hfit)' * l4p;

% Transform points to fit in the new image:
x1fit = Hfit * (H * x1);
x2fit = Hfit * (H * x2);
x3fit = Hfit * (H * x3);
x4fit = Hfit * (H * x4);

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

% Window lines:
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'r');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'g');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'b');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');


% Transformed image:
figure()
imshow(I2)
hold on

% Window corners:
plot(x1fit(1)/x1fit(3), x1fit(2)/x1fit(3), '*r', 'LineWidth', 1)
plot(x2fit(1)/x2fit(3), x2fit(2)/x2fit(3), '*g', 'LineWidth', 1)
plot(x3fit(1)/x3fit(3), x3fit(2)/x3fit(3), '*b', 'LineWidth', 1)
plot(x4fit(1)/x4fit(3), x4fit(2)/x4fit(3), '*y', 'LineWidth', 1)

% Window lines:
t=1:0.1:1000;
plot(t, -(l1fit(1)*t + l1fit(3)) / l1fit(2), 'r');
plot(t, -(l2fit(1)*t + l2fit(3)) / l2fit(2), 'g');
plot(t, -(l3fit(1)*t + l3fit(3)) / l3fit(2), 'b');
plot(t, -(l4fit(1)*t + l4fit(3)) / l4fit(2), 'y');

