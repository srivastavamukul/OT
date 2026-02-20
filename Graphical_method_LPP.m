clc; clear; close all;

%% 1. Problem Setup (Standard Form: Ax <= b)
% Maximize Z = 40x + 30y
% Constraints:
% 1) 2x + 1y <= 10
% 2) 1x + 2y <= 8
% 3) x >= 0, y >= 0 (Handled as lines later)

C = [40, 30];        % Objective Function
A = [2, 1;           % Constraint Matrix
     1, 2];
b = [10; 8];         % Constraint Limits

% Add Non-Negativity explicitly as lines:
% x >= 0 is equivalent to -1x + 0y <= 0
% y >= 0 is equivalent to 0x + -1y <= 0
A_total = [A; -1, 0; 0, -1];
b_total = [b; 0; 0];

%% 2. Find All Intersections (The "Corner" Logic)
num_lines = size(A_total, 1);
points = []; % Store potential corners here

% Loop through every unique pair of lines
for i = 1:num_lines-1
    for j = i+1:num_lines
        % Get coefficients for Line i and Line j
        A_sub = [A_total(i, :); A_total(j, :)];
        b_sub = [b_total(i); b_total(j)];
        
        % Check for parallel lines (Determinant = 0)
        if abs(det(A_sub)) < 1e-6
            continue; % Skip parallel lines
        end
        
        % Solve linear system: A_sub * X = b_sub
        X = A_sub \ b_sub; 
        
        points = [points; X']; % Append to list
    end
end

%% 3. Feasibility Check (The Filter)
feasible_points = [];
tolerance = 1e-6; % To handle floating point errors

for k = 1:size(points, 1)
    pt = points(k, :);
    x_val = pt(1); 
    y_val = pt(2);
    
    % Check against ALL original constraints
    % (A_total * pt') must be <= b_total
    check = (A_total * pt' <= b_total + tolerance);
    
    if all(check)
        feasible_points = [feasible_points; pt];
    end
end

% Remove duplicate points (common in corner intersections)
feasible_points = unique(feasible_points, 'rows');

%% 4. Optimization (The Selection)
if isempty(feasible_points)
    error('No feasible region found!');
end

% Evaluate Z for all feasible points
z_values = (C(1) * feasible_points(:, 1)) + (C(2) * feasible_points(:, 2));

[max_z, idx] = max(z_values);
optimal_pt = feasible_points(idx, :);

%% 5. Output & Plotting (The Visualization)

% A. Print Results
fprintf('---------------------------------\n');
fprintf('Calculated Corners (Feasible):\n');
disp(feasible_points);
fprintf('---------------------------------\n');
fprintf('Optimal Solution Found:\n');
fprintf('x = %.2f\n', optimal_pt(1));
fprintf('y = %.2f\n', optimal_pt(2));
fprintf('Max Profit Z = %.2f\n', max_z);

% B. Plotting Logic
figure; hold on; grid on;
axis([-1 10 -1 10]); % Set axis limits manually for better view
xlabel('X Axis'); ylabel('Y Axis');
title('Graphical Method Implementation (From Scratch)');

% 1. Draw the Feasible Region
% We use 'convhull' to order the points to draw a polygon
k = convhull(feasible_points(:,1), feasible_points(:,2));
fill(feasible_points(k,1), feasible_points(k,2), 'g', 'FaceAlpha', 0.3);

% 2. Draw the Lines (Visualizing boundaries)
x_range = -1:0.1:10;
colors = ['r', 'b', 'k', 'm']; % Colors for different lines

for i = 1:size(A, 1) % Only plot the main structural constraints
    % ax + by = c  =>  by = c - ax  =>  y = (c - ax) / b
    % Handle case where b is 0 (vertical line)
    if A(i, 2) == 0
        x_line = repmat(b(i)/A(i,1), size(x_range));
        plot(x_line, x_range, 'LineWidth', 2, 'DisplayName', sprintf('Constraint %d', i));
    else
        y_line = (b(i) - A(i,1)*x_range) / A(i,2);
        plot(x_range, y_line, 'LineWidth', 2, 'DisplayName', sprintf('Constraint %d', i));
    end
end

% 3. Highlight Optimal Point
plot(optimal_pt(1), optimal_pt(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(optimal_pt(1)+0.2, optimal_pt(2), 'Optimal', 'FontSize', 12, 'FontWeight', 'bold');

legend('Feasible Region');
hold off;
