clc;
clear all;
close all;

hold on;

ang = [0 pi/4]; % Angles in radians, slightly different angles
r = [2 7]; % Vector magnitudes, same size for demonstration
start = [0 0];
resultant = [0 0];

for i = 1:numel(r)
    % Plot the vector
    plot([start(1) start(1) + r(i)*cos(ang(i))], [start(2) start(2) + r(i)*sin(ang(i))], 'b-');
    
    % Plot the arrow
    quiver(start(1), start(2), r(i)*cos(ang(i)), r(i)*sin(ang(i)), 0, 'b');
    
    % Display the magnitude and phase on the vector
    text(start(1) + 0.5*r(i)*cos(ang(i)), start(2) + 0.5*r(i)*sin(ang(i)), ['Magnitude: ' num2str(r(i)) ', \theta: ' num2str(rad2deg(ang(i))) '°'], 'Color', 'b');
    
    % Update the resultant vector
    resultant = resultant + [r(i)*cos(ang(i)) r(i)*sin(ang(i))];
    
    % Update the start point
    start = start + [r(i)*cos(ang(i)) r(i)*sin(ang(i))];
end

% Plot the resultant vector
plot([0 resultant(1)], [0 resultant(2)], 'r-');
quiver(0, 0, resultant(1), resultant(2), 0, 'r');

% Calculate and display the phase angle of the resultant vector
phase_angle = atan2(resultant(2), resultant(1));
text(resultant(1)*0.5, resultant(2)*0.5, ['Magnitude: ' num2str(norm(resultant)) ', \theta: ' num2str(rad2deg(phase_angle)) '°'], 'Color', 'r');

hold off;

% Customize plot appearance
xlabel('Real Axis');
ylabel('Imaginary Axis');
title('Phasor Diagram');

% Adjust axis limits and add margin
margin = 0.5;
xlim([min([0, resultant(1)]) - margin, max([0, resultant(1)]) + margin]);
ylim([min([0, resultant(2)]) - margin, max([0, resultant(2)]) + margin]);

grid on;
legend('Phasor Vectors', 'Resultant Vector');
