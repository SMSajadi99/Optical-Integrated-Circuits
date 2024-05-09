clc; clear all; close all;

x = linspace(-10,10,10^5);
landa = 10;
k = 2*pi / landa;

phi0 = pi / 2;
y0 = sin(k * x + phi0);

phi1 = pi / 3;
y1 = sin(k * x + phi1);

phi2 = 0;
y2 = sin(k * x + phi2);

phi3 = -pi / 2;
y3 = sin(k * x + phi3);

phi4 = 3 * pi / 10;
y4 = sin(k * x + phi4);

figure;
plot(x,y0)
hold on;
grid on;
plot(x,y1)
plot(x,y2)
plot(x,y3)
plot(x,y4)
hold off

xlabel('X (m)')
ylabel('I')
title('y(x) = sin(k * x + phi) ')
legend('y(x) = sin(k * x + pi/2) ','y(x) = sin(k * x + pi/3) ','y(x) = sin(k * x + 0) ','y(x) = sin(k * x - pi/2) ','y(x) = sin(k * x + 3*pi/10) ')

%%

clc; clear all; close all;

x = linspace(-10,10,10^5);
landa = 10;
k = 2*pi / landa;

phi0 = pi / 2 + pi / 2;
y0 = cos(k * x + phi0);

phi1 = pi / 3 + pi / 2;
y1 = cos(k * x + phi1);

phi2 = 0 + pi / 2;
y2 = cos(k * x + phi2);

phi3 = -pi / 2 + pi / 2;
y3 = cos(k * x + phi3);

phi4 = 3 * pi / 10 + pi / 2;
y4 = cos(k * x + phi4);

figure;
plot(x,y0)
hold on;
grid on;
plot(x,y1)
plot(x,y2)
plot(x,y3)
plot(x,y4)
hold off

xlabel('X (m)')
ylabel('I')
title('y(x) = cos(k * x + phi) ')
legend('y(x) = cos(k * x + pi) ','y(x) = cos(k * x + 5*pi/3) ','y(x) = cos(k * x + pi/2) ','y(x) = cos(k * x - 0) ','y(x) = cos(k * x + 8*pi/10) ')

