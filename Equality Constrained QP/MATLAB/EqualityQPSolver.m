%{
    EqualityQPSolver:
    
    Solves Quadratic Program

        min 1/2x'Hx + gx subject to A'x = b, x ≥ 0

    Require:
    - H, g, A, b
    - Solver: Solution method to use
        - 'LDLdense': LDL factorisation with dense KKT matrix
        - 'LDLsparse': LDL factorisation with sparse KKT matrix
        - 'LUdense': LU factorisation with dense KKT matrix
        - 'LUsparse': LU factorisation with sparse KKT matrix
        - 'nullspace': Null-space method
        - 'rangespace': Range-space method
    Ensure:
    - Optimal solution x^* and λ^* found

    Raises error if solver is not recognised.
%}
function [x, lambda] = EqualityQPSolver(H, g, A, b, solver)
    switch solver
        case 'LDLdense'
            [x, lambda] = EqualityQPSolverLDLdense(H, g, A, b);
        case 'LDLsparse'
            [x, lambda] = EqualityQPSolverLDLsparse(H, g, A, b);
        case 'LUdense'
            [x, lambda] = EqualityQPSolverLUdense(H, g, A, b);
        case 'LUsparse'
            [x, lambda] = EqualityQPSolverLUsparse(H, g, A, b);
        case 'nullspace'
            [x, lambda] = EqualityQPSolverNullSpace(H, g, A, b);
        case 'rangespace'
            [x, lambda] = EqualityQPSolverRangeSpace(H, g, A, b);
        otherwise
            error('Solver not recognised');
    end
end
