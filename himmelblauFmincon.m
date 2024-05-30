%{

Solves Himmelblaus test problem, using Matlab's fmincon function.

min f(x)=(x1^2 +x2 -11)^2 +(x1 +x2^2 -7)^2
s.t:
(x1+2)^2 -x2 >= 0, nonlinear inequality constraint
-4x1 +10x2 >= 0, linear inequality constraint
and 
-5 <= x1 <= 5
-5 <= x2 <= 5
%}

% Selection of initial guess
x1 = [-4; 1]; % Good
x2 = [1; 1]; % Good



x0 = x1;

algorithm = 'sqp';
fp = '/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem4/data/fmincon';
[xsol,fval,exitflag,output,lambda,grad,hessian, history] = runfmincon_hb(x0, algorithm);

fprintf('The solution is:\n');
disp(xsol);

plotIter(history, false);


% saving the results 
% filename = strcat(fp, '/results_', algorithm, '.mat');
% writematrix(iterates, filename);

    

% Functions

function [c, ceq, gc, gceq] = himmelnonlcon(X)
    x = X(1);
    y = X(2);
    c = -((x+2)^2 - y);
    ceq = [];
    if nargout > 2
        gc = [-2*(x+2); 1];
        gceq = [];
    end
end



function [xsol,fval,exitflag,output,lambda,grad,hessian, history] = runfmincon_hb(x0, algorithm)
    %{
    Solves Himmelblaus test problem, using Matlab's fmincon function.

    Allows the user to specify the algorithm used and initial guess, plotting and returning iterates.

    %}
    % Set up shared variables with outfun
    history.x = [];
    history.fval = [];
     
    % Defining Himmelblau's problem with constraints

    % Linear inequality constraint, Ax <= b (notice that the inequality is reversed compared to the usual notation)

    A = [4,-10]; 
    b = 0; 

    % No linear equality constraints

    Aeq = []; 
    beq = []; 

    % Lower and upper bounds
    lb = [-5, -5];
    ub = [5, 5];

    options = optimoptions('fmincon',...
    'Display', 'iter', 'Algorithm',algorithm,...
    'SpecifyObjectiveGradient', true,...
    'SpecifyConstraintGradient', true,...
    'HessianFcn', @himmelblauHessian); % SQP and active-set do not use Hessian (ignored)

    [xsol,fval,exitflag,output,lambda,grad,hessian] = fmincon(@himmelblau,x0,A,b,Aeq,beq,lb,ub,@himmelnonlcon, options);
     
end