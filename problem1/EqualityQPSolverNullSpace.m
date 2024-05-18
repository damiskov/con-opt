%{
    EqualityQPSolverNullSpace:
    
    Solves Quadratic Program with dense KKT matrix using Null Space method

        min 1/2x'Hx + gx subject to A'x = b, x ≥ 0

    Require:
    - H, g, A, b 
    Ensure:
    - Optimal solution x^* and λ^* found
%}

function [x, lambda] = EqualityQPSolverNullSpace(H,g,A,b)

[Q, Rbar] = qr(A);
m1 = size(Rbar, 2);
Q1 = Q(:, 1:m1);
Q2 = Q(:, m1+1:end);
R = Rbar(1:m1, 1:m1);

xY = R' \ b;
H_q = Q2' * H * Q2;
Lq = chol(H_q, 'lower');
xZ = -Lq' \ (Lq \ (Q2' * H * Q1 * xY + Q2' * g));
x= Q1 * xY + Q2 * xZ;
lambda = R \ (Q1' * (H * x + g));
end

% function [x, lambda] = EqualityQPSolverNullSpace(H,g,A,b) 
%     [n, ~] = size(A);

%     % Null space method - Orthonormal basis

%     [Q, R] = qr(A);

%     [~, mR] = size(R);

%     Q1 = Q(:, 1:mR);
%     Q2 = Q(:, mR+1:n);
%     Rnew = R(1:mR, 1:mR);

%     % Procedure

%     xY = Rnew' \ b;
%     xZ = -(Q2'*H*Q2) \ (Q2'*H*Q1*xY + Q2'*g);


%     x = Q1*xY + Q2*xZ;
%     lambda = Rnew \ Q1'*(H*x + g);


% end