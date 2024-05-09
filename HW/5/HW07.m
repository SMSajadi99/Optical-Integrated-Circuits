clc; clear all; close all;

x = linspace(-30,30,10^5);
y0 = 4 ./ (x.^2 +2);
y2 = 4 ./ ((x+2.5*2).^2 +2);
y5 = 4 ./ ((x+2.5*5).^2 +2);

figure;
plot(x,y0,'r')
hold on;
grid on;
plot(x,y2,'g')
plot(x,y5,'b')
hold off

xlabel('X (m)')
ylabel('I')
title('y(x,t) = 4 / ((x + 2.5*t)^2 +2)')
legend('y(x,0) = 4 / (x^2 +2)','y(x,2) = 4 / ((x+2.5*2)^2 +2)','y(x,5) = 4 / ((x+2.5*5)^2 +2)')