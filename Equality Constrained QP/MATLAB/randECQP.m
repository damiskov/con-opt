%{
    Generates a "random" Equality Constrained Quadratic Program (ECQP) with n variables.

    Require:
    - n: number of variables (in Objective function)
    Ensure:
    KKT equation matrices:
    - H: n x n symmetric positive definite matrix
    - g: n x 1 vector
    - A: m x n matrix
    - b: n x 1 vector
    Solutions:
    - x: n x 1 vector
    - lambda: m x 1 vector
%}
% function [H, g, A, b, x_opt, lambda_opt] = randECQP(n, sparsity)
    
% if n < 3
%     error('n must be at least 3');
% end
% m = round(rand(1)*n); % m = round(β*n)
% max_val = 5;

% x_opt = max_val*(rand(n, 1)-0.5)*2;

% H = full(max_val*sprandsym(n, sparsity, rand(n, 1)+eps));
% A = max_val*(rand(n, m)-0.5)*2;
% lambda_opt = max_val*(rand(m, 1)-0.5)*2;

% b = A'*x_opt;
% g = A*lambda_opt-H*x_opt;


% end


function [H, g, A, b, x, lambda] = randECQP(n, sparsity)

if n < 3
    error('n must be at least 3');
end

m = round(rand(1)*n); % m = round(β*n)
alpha = 100; % α = rand(1)
M = rand(n, n) .* (rand(n, n) < sparsity/100); % 15% sparsity
H = M' * M + alpha*eye(n);

% A = rand(n,m).* (rand(n,m) < 0.15); % 15% sparsity
A = rand(n, m);

x = rand(n,1);
lambda = rand(m,1);

b = A' * x;

g = A*lambda-H*x;

end