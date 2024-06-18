%{
    EqualityQPSolverRangeSpace:
    
    Solves Quadratic Program with dense KKT matrix using Range Space method

        min 1/2x'Hx + gx subject to A'x = b, x ≥ 0

    - Useful when H is well conditioned and easy to invert.



    Require:
    - H, g, A, b 
    Ensure:
    - Optimal solution x^* and λ^* found
%}

function [x, lambda] = EqualityQPSolverRangeSpace(H,g,A,b)

L = chol(H, "lower");
v = L' \  (L\g);
H_A = A' * (L' \ (L \ A));
L_A = chol(H_A, "lower");
lambda = L_A' \ (L_A \ (A' * v + b));
x = L' \ (L \ (A * lambda - g));

end
