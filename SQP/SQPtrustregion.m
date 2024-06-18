%{
    Performs SQP on general NLP with both equality and inequality constraints. (and bounds)
    - l1 penalty method for QP subproblem.
    - Trust region or line search method for step size.
    - Hessian is updated via BFGS method or given as input.

    Args:
        objective (function): Objective function to be minimized.
        x0 (array): Initial guess.
        x_l (array): Lower bounds for x.
        x_u (array): Upper bounds for x.
        g_l (array): Lower bounds for inequality constraints.
        g_u (array): Upper bounds for inequality constraints.
        h (function): Equality constraint function.
        g (function): Inequality constraint function.
        option (str): Either "trust-region" or "line-search".
         ---- (optional arguments))-----
        maxIter (int): Maximum number of iterations.
        tol (float): Tolerance for convergence.
        objectiveHessian (function): Function to compute Hessian of objective function.
        
        

    Returns:
        x (array): iterate locations.
        f (float): Objective function value at optimal solution.
        exitflag (int): 1 if converged, 0 otherwise.
        iterations (int): Number of iterations performed.
        lambda (array): Lagrange multipliers.
        grad (array): Gradient of objective function at optimal solution.
        H (array): Hessian of objective function at optimal solution.
 
%}

function [x_iter, f_iter, exitflag, iterations, lambda, DxL_norm, H, KKT] = SQPtrustregion(objective, x0, x_l, x_u, g_l, g_u, h, g, objectiveHessian)
    maxIter = 50;
    tol = 1e-6;
    n = length(x0);
    nEq = length(h(x0));
    nIneq = 2*length(g(x0)); % 2x to account for both lower and upper bounds
    k = 0;
    x = zeros(n, maxIter);
    f = zeros(1, maxIter);
    df = zeros(n, maxIter);
    lambdaEq = zeros(nEq, maxIter);
    lambdaIneq = zeros(nIneq, maxIter);
    lambdaLower = zeros(n, maxIter);
    lambdaUpper = zeros(n, maxIter);
    KKT = zeros(1, maxIter);
    DxL_norm = zeros(1, maxIter);
    % Set initial lagrange multipliers to 1
    lambdaEq(:,1) = 5*ones(nEq, 1);
    lambdaIneq(:,1) = 5*ones(nIneq, 1);
    lambdaLower(:,1) = 5*ones(n, 1);
    lambdaUpper(:,1) = 5*ones(n, 1);
    % Check if x0 is within bounds - if not place in the middle of the box
    if any(x0 < x_l) || any(x0 > x_u)
        x0 = (x_l + x_u)/2;
    end
    % Initial values
    x(:,1) = x0;
    [f(1), df(:,1)] = objective(x0);
    [h_k, dh_k, d2h_k] = h(x0);
    [g_k, dg_k, d2G_k] = g(x0); % Writing inequality constraints in standard form
    G_k = [g_k - g_l; g_u - g_k];
    dG_k = [dg_k, -dg_k];
    % Initial Hessian
    H = zeros(n, n, maxIter);
    if nargin < 12
        H(:,:,1) = eye(n); % Default to identity matrix
        hessian_method='BFGS';
    else
        H(:,:,1) = lagrangeHessian(objectiveHessian(x0), lambdaIneq(:,1), lambdaEq(:,1), d2G_k, d2h_k);
        hessian_method='given';
        % Enforcing convexity
        H(:,:,k+1) = regulariseHessian(H(:,:,k+1));
    end

    % Initial Line-Search, trust-region and l1 QP penalty Parameters
    mu = 50;
    % Initial Trust region radius
    Delta_k = 1;
    converged = false;

    % Main loop

    while (k<maxIter && ~converged)

        k = k+1;
        % Box constraints for QP subproblem
        p_k_lb = max(-Delta_k*ones(n, 1), x_l-x(:, k));
        p_k_ub = min(Delta_k*ones(n, 1), x_u-x(:, k));
         
        n_l1 = n + 2*nEq+nIneq;
        H_qp = zeros(n_l1);
        H_qp(1:n, 1:n) = H(:,:,k);
        df_qp = [df(:, k); mu*ones(nEq, 1); mu*ones(nEq, 1); mu*ones(nIneq, 1)];
        A_ineq_qp = [-dG_k', zeros(nIneq, nEq), zeros(nIneq, nEq), -eye(nIneq)];
        b_ineq_qp = G_k;
        A_eq_qp = [dh_k', -eye(nEq), eye(nEq), zeros(nEq, nIneq)];
        b_eq_qp = -h_k;
        x_qp_lb = [p_k_lb; zeros(nEq, 1); zeros(nEq, 1); zeros(nIneq, 1)];
        x_qp_ub = [p_k_ub; 1e6*ones(nEq, 1); 1e6*ones(nEq, 1); 1e6*ones(nIneq, 1)];


        options = optimoptions('quadprog', 'Display', 'off');
        [x_qp, ~, ~, ~, lambda] = quadprog(H_qp, df_qp, A_ineq_qp, b_ineq_qp, A_eq_qp, b_eq_qp, x_qp_lb, x_qp_ub, [], options);

        % Extracting values from l1 QP results
        p_k = x_qp(1:n);
        lambda.ineqlin = lambda.ineqlin(1:nIneq);
        lambda.eqlin = lambda.eqlin(1:nEq);
        lambda.lower = lambda.lower(1:n);
        lambda.upper = lambda.upper(1:n);

        p_lambdaEq = lambda.eqlin-lambdaEq(:, k);
        p_lambdaIneq = lambda.ineqlin-lambdaIneq(:, k);
        p_lambdaLower = lambda.lower-lambdaLower(:, k);
        p_lambdaUpper =  lambda.upper-lambdaUpper(:, k);

    
    
        alpha = 1;
        

        x_test = x(:, k) + alpha*p_k;
        [f_test, df_test] = objective(x_test);
        [h_test, dh_test, ~] = h(x_test);
        [g_test, dg_test, ~] = g(x_test);
        G_test = [g_test - g_l; g_u - g_test];
        dG_test = [dg_test, -dg_test];

        
        if norm(p_k, 'inf') < tol
            rho = 1;
        else   
            rho = (f_test - f(k))/(df(:, k)'*p_k + 1/2*p_k'*H(:,:,k)*p_k);
        end

        if rho < 0.25
            Delta_k = 0.25*Delta_k;
        elseif rho > 0.25 && rho <= 0.75
            Delta_k = min(2*Delta_k, 1);
        end

        if rho > 0
            take_step = true;
        else
            take_step = false;
        end
        


        if take_step % The Step has been accepted by trust-region

            % Taking step
            x(:, k+1) = x(:, k) + alpha*p_k;
            lambdaEq(:, k+1) = lambdaEq(:, k) + alpha*p_lambdaEq;
            lambdaIneq(:, k+1) = lambdaIneq(:, k) + alpha*p_lambdaIneq;
            lambdaLower(:, k+1) = lambdaLower(:, k) + alpha*p_lambdaLower;
            lambdaUpper(:, k+1) = lambdaUpper(:, k) + alpha*p_lambdaUpper;

            D_x_Lagrange_1 = df(:, k) -  dh_k*lambdaEq(:, k+1) - dG_k*lambdaIneq(:, k+1) - lambdaLower(:, k+1) + lambdaUpper(:, k+1);
     
            [f(k+1), df(:, k+1)] = objective(x(:, k+1));
            [h_k, dh_k, d2h_k] = h(x(:, k+1));
            [g_k, dg_k, d2G_k] = g(x(:, k+1)); % Writing inequality constraints in standard form
            G_k = [g_k - g_l; g_u - g_k];
            dG_k = [dg_k, -dg_k];

            D_x_Lagrange_2 = df(:, k+1) -  dh_k*lambdaEq(:, k+1) - dG_k*lambdaIneq(:, k+1) - lambdaLower(:, k+1) + lambdaUpper(:, k+1);

            p = x(:, k+1) - x(:, k); % Should be equal to alpha*p_k
            q = D_x_Lagrange_2 - D_x_Lagrange_1;
            DxL_norm(:, k) = norm(D_x_Lagrange_2, 'inf');
            KKT(:, k) = norm(D_x_Lagrange_2, 'inf') + norm(G_k.*lambdaIneq(:, k+1), 'inf') + norm(h_k.*lambdaEq(:, k+1), 'inf') + norm((x_l-x(:, k+1)).*lambdaLower(:, k+1), 'inf') + norm((x_u-x(:, k+1)).*lambdaUpper(:, k+1), 'inf');

            % Update Hessian

            if strcmp(hessian_method, 'BFGS')
                % BFGS update
                B = H(:,:,k);
                theta = 1;
                if p'*q < 0.2*p'*B*p
                    theta = 0.8*p'*(B*p)/(p'*(B*p) - p'*q);
                end
                
                r = theta*q + (1-theta)*B*p;

                Bk = B - (B*p)*(p'*B)/(p'*B*p) + (r*r')/(p'*r);
                H(:,:,k+1) = Bk;
                % Enforcing convexity
                H(:,:,k+1) = regulariseHessian(H(:,:,k+1));
                
            else
                % Using given Objective Hessian
                H(:,:,k+1) = lagrangeHessian(objectiveHessian(x(:, k+1)), lambdaIneq(:, k+1), lambdaEq(:, k+1), d2G_k, d2h_k);
                % Enforcing convexity
                H(:,:,k+1) = regulariseHessian(H(:,:,k+1));
                

            end

            % Check convergence

            converged = norm(D_x_Lagrange_2, 'inf') <= tol &&...
            norm(G_k.*lambdaIneq(:, k+1), 'inf') <= tol &&...
            norm(h_k.*lambdaEq(:, k+1), 'inf') <= tol &&...
            norm((x_l-x(:, k+1)).*lambdaLower(:, k+1), 'inf') <= tol &&...
            norm((x_u-x(:, k+1)).*lambdaUpper(:, k+1), 'inf') <= tol;
            fprintf('Equality constraint lagrange multiplier:\n');
            disp(lambdaEq(:, k+1));

            % Check if Hessian contains Nan
            if any(isnan(H(:,:,k+1)))
                fprintf('Hessian contains NaN\n');
                break;
            end

        else
            %{
                This section is only reached when trust-region radius must be increased.
                - Still counts as an iteration.
                - Values are not updated.
            %}
            x(:, k+1) = x(:, k);
            lambdaEq(:, k+1) = lambdaEq(:, k);
            lambdaIneq(:, k+1) = lambdaIneq(:, k);
            lambdaLower(:, k+1) = lambdaLower(:, k);
            lambdaUpper(:, k+1) = lambdaUpper(:, k);
            f(k+1) = f(k);
            df(:, k+1) = df(:, k);
            H(:,:,k+1) = H(:,:,k);

        end
    end

    % Converged or reached maximum number of iterations
    if converged
        exitflag = 1;
        fprintf('Converged after %d iterations\n', k);
    else
        exitflag = 0;
        fprintf('Did not converge after %d iterations\nEncountered error.', k);
    end
    iterations = k;
    x_iter = x(:, 1:k)';
    lambda.eq = lambdaEq(:, k);
    lambda.ineq = lambdaIneq(:, k);
    lambda.lower = lambdaLower(:, k);
    lambda.upper = lambdaUpper(:, k);
    DxL_norm = DxL_norm(:, 1:k)';
    KKT = KKT(:, 1:k)';
    f_iter = f(1:k)';
    grad = df;
    % Output results of the SQP algorithm
    % Output results of the SQP algorithm
    fprintf('-------- Trust-Region SQP Complete ---------\n');
    fprintf('Hessian method: %s\n', hessian_method);
    fprintf('Initial Point: \n');
    disp(x0);
    if exitflag == 1
        fprintf('Program successfully converged\n');
    else
        fprintf('Program did not converge\n');
    end
    fprintf('\nOptimal solution: \n');
    disp(x(:,k));
    fprintf('Optimal function value: \n');
    disp(f(k));
end


function Lxx = lagrangeHessian(H, lambdaIneq, lambdaEq, d2G, d2h)
    %{
        Calculates the Hessian of the Lagrangian function.
        
        Used for QP subproblem.
    %}

    sumh = 0;
    sumG = 0;
    for i = 1:length(lambdaEq)
        sumh = sumh + lambdaEq(i)*d2h(:, i);
    end
    for i = 1:length(lambdaIneq)
        sumG = sumG + lambdaIneq(i)*d2G(:, i);
    end

    Lxx = H - sumG - sumh;

end

function H_reg = regulariseHessian(H)

    eigvals = eig(H);
    lambda_min = min(eigvals);

    if lambda_min <= 0
        alpha = -lambda_min + 1e-6;
    else
        alpha = 0;
    end
    
    % Modify the Hessian
    H_reg = H + alpha * eye(size(H));
end