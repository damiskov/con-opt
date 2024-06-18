%{
    EqualityQPSolverLUsparse:
    
    Solves Quadratic Program with sparse KKT matrix using slightly modified LU factorisation

        min 1/2x'Hx + gx subject to A'x = b, x ≥ 0

    Require:
    - H, g, A, b 
    Ensure:
    - Optimal solution x^* and λ^* found
%}

function [x, lambda] = EqualityQPSolverLUsparse(H,g,A,b)
    
[n, m] = size(A);
KKT = [[H, -A]; [-A', zeros(m)]];
v = - [g; b];

% Make KKT and RHS sparse
KKT = sparse(KKT);
v = sparse(v);


[L, U, P] = lu(KKT, "vector");
sol = U \ (L \ (v(P)));

lambda = sol(n+1:end);
x = sol(1:n); 
    
end