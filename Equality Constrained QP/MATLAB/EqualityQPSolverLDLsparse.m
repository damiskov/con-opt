%{
    EqualityQPSolverLDLdense:
    
    Solves Quadratic Program with dense KKT matrix using LDL factorisation

        min 1/2x'Hx + gx subject to A'x = b, x ≥ 0

    Require:
    - H, g, A, b 
    Ensure:
    - Optimal solution x^* and λ^* found
%}

function [x, lambda] = EqualityQPSolverLDLsparse(H,g,A,b)
[n, m] = size(A);
KKT = [[H, -A]; [-A', zeros(m)]];
v = - [g; b];

% Convert KKT to sparse matrix
KKT = sparse(KKT);
v = sparse(v);

% LDL factorisation

[L, D, P] = ldl(KKT);
x = P*(L'\(D\(L\(P'*v))));

lambda = x(n+1:end);
x = x(1:n); 
end