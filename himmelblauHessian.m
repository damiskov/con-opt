function H = himmelblauHessian(x)
%{
    Hessian of the objective function for our Himmelblau test problem.
%}

H = [12*x(1)^2 + 4*x(2) - 42, 4*x(1) + 4*x(2);
    4*x(1) + 4*x(2), 4*x(1) + 12*x(2)^2 - 26];
end