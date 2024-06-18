H = [5.0000 1.8600 1.2400 1.4800 -0.4600;
     1.8600 3.0000 0.4400 1.1200 0.5200;
     1.2400 0.4400 3.8000 1.5600 -0.5400;
     1.4800 1.1200 1.5600 7.2000 -1.1200;
    -0.4600 0.5200 -0.5400 -1.1200 7.8000];

g = [-16.1000;
    -8.5000;
    -15.7000;
    -10.0200;
    -18.6800];


A = [16.1000 1.0000; 
     8.5000 1.0000; 
     15.7000 1.0000; 
     10.0200 1.0000; 
     18.6800 1.0000];

% b = [15;
%     1];

% Altering b(1) from 8.5 to 18.68 with n values in between
n = 1000;
b1_vals = linspace(8.5, 18.68, n);

b = [linspace(8.5, 18.68, n);
    ones(1, n)];

% storing solutions for each solver
x_LDLdense = zeros(5,n);
lambda_LDLdense = zeros(2,n);

x_LDLsparse = zeros(5,n);
lambda_LDLsparse = zeros(2,n);

x_LUdense = zeros(5,n);
lambda_LUdense = zeros(2,n);

x_LUsparse = zeros(5,n);
lambda_LUsparse = zeros(2,n);

x_nullspace = zeros(5,n);
lambda_nullspace = zeros(2,n);

x_rangespace = zeros(5,n);
lambda_rangespace = zeros(2,n);

obj_vals = zeros(1, n);


% Solve using dense LDL
for i = 1:n
    [x_LDLdense(:,i), lambda_LDLdense(:,i)] = EqualityQPSolver(H, g, A, b(:,i), 'LDLdense');
    [x_LDLsparse(:,i), lambda_LDLsparse(:,i)] = EqualityQPSolver(H, g, A, b(:,i), 'LDLsparse');
    [x_LUdense(:,i), lambda_LUdense(:,i)] = EqualityQPSolver(H, g, A, b(:,i), 'LUdense');
    [x_LUsparse(:,i), lambda_LUsparse(:,i)] = EqualityQPSolver(H, g, A, b(:,i), 'LUsparse');
    [x_nullspace(:,i), lambda_nullspace(:,i)] = EqualityQPSolver(H, g, A, b(:,i), 'nullspace');
    [x_rangespace(:,i), lambda_rangespace(:,i)] = EqualityQPSolver(H, g, A, b(:,i), 'rangespace');


    % Using LU factorisation solution for objective function value (seems to have the lowest error)
    x = x_LUdense(:,i);
    obj_vals(:,i) = objective_func(x, H, g);

end


% Saving data

save('/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem1/data/x_LDLdense.mat', 'x_LDLdense');
save('/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem1/data/lambda_LDLdense.mat', 'lambda_LDLdense');
save('/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem1/data/x_LDLsparse.mat', 'x_LDLsparse');
save('/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem1/data/lambda_LDLsparse.mat', 'lambda_LDLsparse');
save('/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem1/data/x_LUdense.mat', 'x_LUdense');
save('/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem1/data/lambda_LUdense.mat', 'lambda_LUdense');
save('/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem1/data/x_LUsparse.mat', 'x_LUsparse');
save('/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem1/data/lambda_LUsparse.mat', 'lambda_LUsparse');
save('/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem1/data/x_nullspace.mat', 'x_nullspace');
save('/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem1/data/lambda_nullspace.mat', 'lambda_nullspace');
save('/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem1/data/x_rangespace.mat', 'x_rangespace');
save('/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem1/data/lambda_rangespace.mat', 'lambda_rangespace');

% Saving objective values
save('/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem1/data/b1_vals.mat', 'b1_vals');
save('/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem1/data/obj_vals.mat', 'obj_vals');




function val = objective_func(x, H, g)
    val = 1/2 * x' * H * x + g'*x;
end 