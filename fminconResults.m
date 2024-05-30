%{

Solves Himmelblaus test problem, using Matlab's fmincon function.

min f(x)=(x1^2 +x2 -11)^2 +(x1 +x2^2 -7)^2
s.t:
(x1+2)^2 -x2 >= 0, nonlinear inequality constraint
-4x1 +10x2 >= 0, linear inequality constraint
and 
-5 <= x1 <= 5
-5 <= x2 <= 5

Top left minima:
- x0 = [-5, 1.4]

Bottom left minima:
x0 = [-5, 0]

Global:
x0 = [1, 0]
%}

% Selection of current starting points

x1 =  [-2; 0.8];
x2 = [-5; 0];
x3 =  [1; 0];
x4 = [1; 5];


algorithm = 'sqp';
fp = '/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem4/SQPfinal/fmincon/';

start_points = {x1, x2, x3, x4};
for i = 1:length(start_points)
    x0 = start_points{i};
    disp(x0)
    [xsol,fval,exitflag,output,lambda,grad,hessian, iterates] = runfmincon_hb(x0, algorithm);
    disp(iterates)
    plotIter(iterates', false);
    save_figure(gcf, [fp, 'fmincon_', num2str(i)], 300, 6, 4);
    save([fp, 'fmincon_', num2str(i), '.mat'], 'xsol', 'fval', 'exitflag', 'output', 'lambda', 'grad', 'hessian', 'iterates');
end


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



function [xsol,fval,exitflag,output,lambda,grad,hessian, iterates] = runfmincon_hb(x0, algorithm)
    %{
    Solves Himmelblaus test problem, using Matlab's fmincon function.

    %}
    
    iterates = []; % Set up shared variables with outfun

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
    'OutputFcn', @outputFcn,...
    'SpecifyObjectiveGradient', true,...
    'SpecifyConstraintGradient', true,...
    'HessianFcn', @himmelblauHessian); % SQP and active-set do not use Hessian (ignored)

    [xsol,fval,exitflag,output,lambda,grad,hessian] = fmincon(@himmelblau,x0,A,b,Aeq,beq,lb,ub,@himmelnonlcon, options);

    iterates = evalin('base', 'iterates');
     
end

function stop = outputFcn(x, optimValues, state)
    persistent iterates;
    stop = false;

    switch state
        case 'init'
            iterates = x';
        case 'iter'
            iterates = [iterates; x'];
        case 'done'
            % Save iterates to a .mat file or assign to a variable in the base workspace
            assignin('base', 'iterates', iterates);
    end
end
