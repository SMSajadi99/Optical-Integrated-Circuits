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

[px,py] = gradient(z);

figure
contour(x,y,z)
hold on
quiver(x,y,px,py)
hold off

colorbar

%%
% Define the scalar function
f = @(x, y) sin(x) .* cos(2*y);

% Define the grid
x = linspace(-pi, pi, 100);
y = linspace(-pi, pi, 100);
[X, Y] = meshgrid(x, y);

% Compute the gradient of the scalar field
[grad_x, grad_y] = gradient(f(X, Y), x, y);

% Plot the vector field and its curl
figure;

% Plot the vector field
quiver(X, Y, grad_x, grad_y, 'AutoScale', 'on');
title('Gradient of Scalar Field');

