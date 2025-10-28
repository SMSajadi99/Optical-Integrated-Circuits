clc; clear all; close all;

% Parameters
r = 3.1e-6; % in meters (converted from micrometers)
kappa = 0.52;
a = 0.83;
lambda = linspace(1540e-9, 1560e-9, 1e3); % wavelength range from 500 nm to 600 nm
n_eff = 3.48; % effective refractive index

% Theta calculation function
theta = @(n_eff, r, lambda) 4 * pi^2 * n_eff^2 * (r ./ lambda);

% Define phi and t
phi = 0;
t = sqrt(1 - kappa^2);

% Function to plot
P_B = @(t, theta, phi, a) (a^2 + t^2 - 2*a*t*cos(theta + phi)) ./ (1 + a^2*t^2 - 2*a*t*cos(theta + phi));

P_D = @(t, theta, phi, a) (a^2 * (1 - t^2)) ./ (1 + a^2*t^2 - 2*a*t*cos(theta + phi));


% Calculate theta
theta_val = theta(n_eff, r, lambda);

% Plot the function
figure;
hold on;
grid on;
plot(lambda, P_B(t, theta_val, phi, a));
plot(lambda, P_D(t, theta_val, phi, a));
xlabel('Wavelength (m)');
ylabel('Indensity (a.u.)');
title('Plot of the function Indensity over Wavelength');
legend("P_B","P_D")