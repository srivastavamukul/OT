clc;
clear;

% Min Z = 2x1 + 3x2
% x1 + x2 >= 4
% x1 + 2x2 = 6
% x1, x2 >= 0

M = 1000;

% Tableau: x1 x2 s1 a1 a2 | RHS
% Bottom row stores reduced costs (c_j - z_j), initialized to +c_j.
% Surplus s1 has cost 0; artificials a1, a2 have cost +M.
T = [1 1 -1 1 0 4;
     1 2  0 0 1 6;
     2 3  0 M M 0];

[m, n] = size(T);

% Drive artificials out of the objective row.
% a1 is basic in row 1, a2 is basic in row 2. Subtract M*row to zero them.
T(end,:) = T(end,:) - M * T(1,:);
T(end,:) = T(end,:) - M * T(2,:);

max_iter = 50;
iter = 0;

while true
    iter = iter + 1;
    if iter > max_iter
        disp('Iteration limit reached'); break;
    end

    % For MIN: entering variable = most negative reduced cost
    [min_val, pivot_col] = min(T(end, 1:end-1));
    if min_val >= -1e-9
        break;
    end

    % Leaving variable: minimum ratio test
    ratios = inf(m-1, 1);
    for i = 1:m-1
        if T(i, pivot_col) > 1e-9
            ratios(i) = T(i, end) / T(i, pivot_col);
        end
    end
    [min_ratio, pivot_row] = min(ratios);
    if isinf(min_ratio)
        disp('Unbounded solution'); break;
    end

    % Pivot
    T(pivot_row, :) = T(pivot_row, :) / T(pivot_row, pivot_col);
    for i = 1:m
        if i ~= pivot_row
            T(i, :) = T(i, :) - T(i, pivot_col) * T(pivot_row, :);
        end
    end
end

disp('Final Tableau:');
disp(T);

disp('Optimal value of Z:');
disp(-T(end, end));   % Z = -bottom-right in this convention
