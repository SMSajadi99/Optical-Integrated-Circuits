% 🌟✨🖥️ MATLAB Code by Seyed Mohammad Sajadi 🖥️✨🌟
% Student Number: 402448040

% Define the function
f = @(x, y) sin(x) .* cos(2*y);

% Generate x and y values
x = linspace(-pi, pi, 100);
y = linspace(-pi, pi, 100);

% Create a meshgrid from x and y
[x, y] = meshgrid(x, y);

% Evaluate the function at each point in the meshgrid
z = f(x, y);

% Create a 3D plot
figure;
surf(x, y, z);
title('😊 Function Plot: f(x, y) = sin(x)*cos(2*y) 😊');
xlabel('x');
ylabel('y');
zlabel('f(x, y)');

colorbar
