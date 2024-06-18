%{
    EqualityQPSolverLUdense:
    
    Solves Quadratic Program with dense KKT matrix using standard LU factorisation and back substitution

        min 1/2x'Hx + gx subject to A'x = b, x ≥ 0

    
    Require:
    - H, g, A, b 
    Ensure:
    - Optimal solution x^* and λ^*  found
%}

function [x, lambda] = EqualityQPSolverLUdense(H,g,A,b)
    
[n, m] = size(A);
KKT = [[H, -A]; [-A', zeros(m)]];
v = - [g; b];


[L, U, P] = lu(KKT, "vector");
sol = U \ (L \ (v(P)));

lambda = sol(n+1:end);
x = sol(1:n); 

end