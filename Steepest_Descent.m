clc;
clear;

% Function
f = @(x1,x2) x1.^2 + 2*x2.^2;

% Gradient
grad_f = @(x1,x2) [2*x1; 4*x2];

% Hessian (constant for this quadratic)
H = [2 0; 0 4];

% Initial point
x = [3; 3];

max_iter = 100;
tol = 1e-6;

fprintf('Initial point: (%f, %f)\n', x(1), x(2));

for iter = 1:max_iter
    
    gradient = grad_f(x(1), x(2));
    
    % Stopping condition
    if norm(gradient) < tol
        fprintf('Converged after %d iterations\n', iter-1);
        break;
    end
    
    % Optimal step size (steepest descent)
    alpha = (gradient' * gradient) / (gradient' * H * gradient);
    
    % Update rule
    x = x - alpha * gradient;
    
    % Print progress
    fprintf('Iter %d: x = (%f, %f), f = %f\n', iter, x(1), x(2), f(x(1), x(2)));
end

fprintf('Optimal point: (%f, %f)\n', x(1), x(2));
fprintf('Minimum value: %f\n', f(x(1), x(2)));
