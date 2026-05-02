clc;
clear;

% Min Z = 8x11 + 6x12 + 10x13 + 9x21 + 7x22 + 4x23
% Supply: S1 = 20, S2 = 30
% Demand: D1 = 10, D2 = 25, D3 = 15

% Cost matrix
C = [8 6 10;
     9 7 4];

% Supply and Demand
supply = [20 30];
demand = [10 25 15];

[m, n] = size(C);
allocation = zeros(m,n);

while any(supply > 0) && any(demand > 0)

    % Step 1: Find minimum cost
    min_cost = inf;

    for i = 1:m
        for j = 1:n
            if supply(i) > 0 && demand(j) > 0
                if C(i,j) < min_cost
                    min_cost = C(i,j);
                    row = i;
                    col = j;
                end
            end
        end
    end

    % Step 2: Allocate
    x = min(supply(row), demand(col));
    allocation(row, col) = x;

    % Step 3: Update
    supply(row) = supply(row) - x;
    demand(col) = demand(col) - x;
end

disp('Allocation Matrix:');
disp(allocation);

total_cost = sum(sum(allocation .* C));

disp('Total Transportation Cost:');
disp(total_cost);
