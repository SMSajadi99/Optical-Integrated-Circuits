% Welcome to the Exciting MATLAB World!
% 🌟 Practice Session: HW - 06 🚀
% 🧠 Your Brain Gym: Tackling MATLAB with Seyed Mohammad Sajadi's Magic
% 🎯 Student Number: 402448040
% Ready, set, code! 💻✨

% Clear the workspace and command window
close all
clc

% Initialize input values (commented out for hardcoded values)
% d = input("Enter the L: ");
% n1 = input("Enter the n1: ");
% n2 = input("Enter the n2: ");
% landa = input("Enter the landa: ");
% alpha = input("Enter the alpha: ");
% theta = input("Enter the theta: ");
n1 = 1.6;
n2 = 1.5;
landa = 530 * 10^-9;
d = 800 * 10^-9;
phi = 4 * pi;

% Calculate the radius of the circle
radius = (2*pi/landa)* d * sqrt(n1^2 - n2^2);

% Display the calculated radius
disp(['The calculated radius of the circle is: ', num2str(radius)]);

% Generate x values and compute y values for tan(x) and cot(x)
x = linspace(0, phi, 10^8);
y_tan =  (n2 / n1)^2 * x .* tan(x);
y_cot = -(n2 / n1)^2 * x .* cot(x);

% Plot the functions
figure;
plot(x, y_tan);
hold on; % Maintain the current plot for additional plots
plot(x, y_cot);
xlim([0 phi]);
ylim([0 10]); % Adjust y-axis limits for better visualization

% Plot the circle using parametric equations
theta = linspace(0, phi, 1000000);
x_circle = radius * cos(theta);
y_circle = radius * sin(theta);
plot(x_circle, y_circle);

% Iterate through x values and find intersection points
intersection_points_tan = [];
intersection_points_cot = [];

tolerance = 0.01; % Tolerance for considering points close

for i = 1:length(x)
    % Check if the point is on the circle
    if abs(x(i)^2 + y_tan(i)^2 - radius^2) < tolerance
        % Check if the point is close to any existing point in the array
        is_close = false;
        for j = 1:size(intersection_points_tan, 1)
            if norm([x(i), y_tan(i)] - intersection_points_tan(j, :)) < tolerance
                is_close = true;
                break;
            end
        end
        
        % If the point is not close to any existing point, add it to the array
        if ~is_close
            intersection_points_tan = [intersection_points_tan; x(i), y_tan(i)];
        end
    end

 
    % Check if the point is on the circle
    if abs(x(i)^2 + y_cot(i)^2 - radius^2) < tolerance
        % Check if the point is close to any existing point in the array
        is_close = false;
        for j = 1:size(intersection_points_cot, 1)
            if norm([x(i), y_cot(i)] - intersection_points_cot(j, :)) < tolerance
                is_close = true;
                break;
            end
        end
        
        % If the point is not close to any existing point, add it to the array
        if ~is_close
            intersection_points_cot = [intersection_points_cot; x(i), y_cot(i)];
        end
    end
end


% Display only positive intersection points with y = x*tan(x)
disp('Intersection points with y = (n2 / n1)^2 * x*tan(x):');
positive_intersection_tan = intersection_points_tan(intersection_points_tan(:, 2) > 0, :);
disp(positive_intersection_tan);

% Display only positive intersection points with y = -x*cot(x)
disp('Intersection points with y = -(n2 / n1)^2 * x*cot(x):');
positive_intersection_cot = intersection_points_cot(intersection_points_cot(:, 2) > 0, :);
disp(positive_intersection_cot);


% Count the number of points for tangent and cotangent
num_points_tan = size(positive_intersection_tan, 1);
num_points_cot = size(positive_intersection_cot, 1);

% Display the counts
disp(['Number of points for tangent: ', num2str(num_points_tan)]);
disp(['Number of points for cotangent: ', num2str(num_points_cot)]);
disp(['Total number of points: ', num2str(num_points_tan + num_points_cot)]);


% Plot intersection points
if ~isempty(intersection_points_tan)
    scatter(intersection_points_tan(:, 1), intersection_points_tan(:, 2), 'r', 'filled');
end

if ~isempty(positive_intersection_cot)
    scatter(positive_intersection_cot(:, 1), positive_intersection_cot(:, 2), 'b', 'filled');
end

% Add labels, title, and legend for clarity
xlabel('x');
ylabel('y');
title('Plot of y = (n2 / n1)^2 * x * tan(x), y = -(n2 / n1)^2 * x * cot(x), and a circle with radius r');
legend('y = (n2 / n1)^2 * x * tan(x)', 'y = -(n2 / n1)^2 * x * cot(x)', 'Circle', 'Intersection Points with tan(x)', 'Positive Intersection Points with cot(x)');
grid on;
hold off;


new_positive_intersection_tan = positive_intersection_tan / d;
new_positive_intersection_cot = positive_intersection_cot / d;

% Define the x values for the entire range
x = linspace(-2 * d, 2 * d, 100000); % Adjust the number of points as needed

% Define the corresponding y values for each function
y = zeros(size(x));

row_positive_intersection_cot = size(positive_intersection_cot, 1);
row_positive_intersection_tan = size(positive_intersection_tan, 1);

num_rows = ceil((row_positive_intersection_tan + row_positive_intersection_cot) / 2);
figure;


for i = 1 : (row_positive_intersection_tan + row_positive_intersection_cot)
    subplot(num_rows, 2, i); % Adjust the layout dynamically
    
    if mod(i, 2) == 0
        % Even index
        K = new_positive_intersection_cot(i / 2, 1);
        k = new_positive_intersection_cot(i / 2, 2);

        for j = 1:length(x)
            % First interval: 0 to d
            if x(j) >= 0 && x(j) <= d
                y(j) = sin(K * x(j));
            end
            
            % Second interval: L to 2*d
            if x(j) > d && x(j) <= 2*d
                if sin(K * d) > 0
                    c = sin(K * d) / exp(-k * d);
                    y(j) = c * exp(-k * x(j));
                else
                    c = -sin(K * d) / exp(-k * d);
                    y(j) = -c * exp(-k * x(j));
                end
            end
            
            % Third interval: 0 to -d
            if x(j) >= -d && x(j) < 0
                y(j) = sin(K * x(j));
            end
            
            % Fourth interval: -d to -2*d
            if x(j) >= -2 * d && x(j) < -d
                if sin(K * (-d)) > 0
                    c = sin(K * -d) / exp(k * -d);
                    y(j) = c * exp(k * x(j));
                else
                    c = -sin(K * -d) / exp(k * -d);
                    y(j) = -c * exp(k * x(j));
                end
            end
        end
        

    else
        % Odd index
        K = new_positive_intersection_tan((i+1) / 2, 1);
        k = new_positive_intersection_tan((i+1) / 2, 2);

        for j = 1:length(x)
            % First interval: 0 to d
            if x(j) >= 0 && x(j) <= d
                y(j) = cos(K * x(j));
            end
            
            % Second interval: L to 2 * d
            if x(j) > d && x(j) <= 2 * d
                if cos(K * d) > 0
                    c = cos(K * d) / exp(-k * d);
                    y(j) = c * exp(-k * x(j));
                else
                    c = -cos(K * d) / exp(-k * d);
                    y(j) = -c * exp(-k * x(j));
                end
            end
            
            % Third interval: 0 to -d
            if x(j) >= -d && x(j) < 0
                y(j) = cos(K * x(j));
            end
            
            % Fourth interval: -d to -2*d
            if x(j) >= -2*d && x(j) < -d
                if cos(K * (-d)) > 0
                    c = cos(K * -d) / exp(k * -d);
                    y(j) = c * exp(k * x(j));
                else
                    c = -cos(K * -d) / exp(k * -d);
                    y(j) = -c * exp(k * x(j));
                end
            end
        end

    end
    

    plot(x, y, 'b-', 'LineWidth', 2);
    
    if mod(i, 2) == 0
        title(['Even Subplot ' num2str(i/2)]);
    else
        title(['Odd Subplot ' num2str((i+1)/2)]);
    end
    
    xlabel('x');
    ylabel('ψ(x)');
end
