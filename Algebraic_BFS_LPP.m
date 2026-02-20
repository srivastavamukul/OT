% --- ALGEBRAIC METHOD LPP (HARD-CODED) ---

% 1. Define the hard-coded matrices for standard form
C = [2, 3, 0, 0];       % Objective function coefficients (including slack)
A = [1, 1, 1, 0;        % Coefficient matrix from constraints
     2, 1, 0, 1];
b = [4;                 % Right-hand side values
     5];

% 2. Identify dimensions
n = size(A, 2);         % Total number of variables (4)
m = size(A, 1);         % Total number of constraints/equations (2)

% 3. Find all combinations of basic variables
% We need to choose 'm' variables out of 'n' to be our basic variables.
% The remaining n-m variables will be set to 0.
combos = nchoosek(1:n, m); 
num_combos = size(combos, 1);

% Initialize variables to keep track of the best solution
best_Z = -inf;
best_solution = zeros(n, 1);

fprintf('--- Calculating Basic Solutions ---\n');

% 4. Loop through every combination
for i = 1:num_combos
    
    % Get the indices of the basic variables for this iteration
    basic_indices = combos(i, :);
    
    % Extract the square matrix (B) for these basic variables
    B = A(:, basic_indices);
    
    % Check if the matrix B is invertible (determinant is not zero)
    % If det(B) == 0, the lines are parallel and there is no intersection here
    if det(B) ~= 0 
        
        % Solve the system of equations B * x_b = b  =>  x_b = B^-1 * b
        % We use MATLAB's \ operator for solving the linear system manually
        x_b = B \ b; 
        
        % Create a full solution vector of zeros
        solution = zeros(n, 1);
        
        % Place the solved basic variable values into their correct spots
        solution(basic_indices) = x_b;
        
        % 5. Feasibility Check: Are all variables >= 0?
        % We use a tiny tolerance (1e-7) to account for floating-point errors
        if all(solution >= -1e-7) 
            
            % 6. Calculate the Objective Function value (Z)
            Z = C * solution;
            
            fprintf('Basic Feasible Solution found at: [%.2f, %.2f, %.2f, %.2f] -> Z = %.2f\n', ...
                    solution(1), solution(2), solution(3), solution(4), Z);
            
            % Check if this is our best (maximum) Z so far
            if Z > best_Z
                best_Z = Z;
                best_solution = solution;
            end
        else
            fprintf('Infeasible Solution (Negative values) at basic indices [%d, %d]\n', basic_indices(1), basic_indices(2));
        end
    else
        fprintf('No intersection (Singular Matrix) at basic indices [%d, %d]\n', basic_indices(1), basic_indices(2));
    end
end

% 7. Print the final optimal solution
fprintf('\n--- OPTIMAL SOLUTION ---\n');
fprintf('Maximum Z = %.2f\n', best_Z);
fprintf('Achieved at x1 = %.2f, x2 = %.2f\n', best_solution(1), best_solution(2));
fprintf('Slack variables: s1 = %.2f, s2 = %.2f\n', best_solution(3), best_solution(4));
