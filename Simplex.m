clc;
clear;

% Max Z = 3x1 + 5x2
% x1 + 2x2 <= 6
% 3x1 + 2x2 <= 12

% Initial Tableau
% [x1 x2 s1 s2 | RHS]

T = [1 2 1 0 6;
     3 2 0 1 12;
    -3 -5 0 0 0];

[m, n] = size(T);

while true
    % Step 1: Entering variable
    [min_val, pivot_col] = min(T(end,1:end-1));

    if min_val >= 0
        break;
    end

    % Step 2: Leaving variable
    ratios = inf(m-1,1);
    for i = 1:m-1
        if T(i,pivot_col) > 0
            ratios(i) = T(i,end) / T(i,pivot_col);
        end
    end

    [~, pivot_row] = min(ratios);

    % Step 3: Pivot
    pivot_element = T(pivot_row, pivot_col);
    T(pivot_row,:) = T(pivot_row,:) / pivot_element;

    for i = 1:m
        if i ~= pivot_row
            T(i,:) = T(i,:) - T(i,pivot_col) * T(pivot_row,:);
        end
    end
end

disp('Final Tableau:');
disp(T);

disp('Optimal value of Z:');
disp(T(end,end));
