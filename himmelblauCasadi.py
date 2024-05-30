import casadi as ca
import numpy as np

x = ca.MX.sym('x', 2)

f = (x[0]**2 + x[1] - 11)**2 + (x[0] + x[1]**2 - 7)**2


# Define inequality constraints
g1 = (x[0] + 2)**2 - x[1]
g2 = -4*x[0] + 10*x[1]

# Add box constraints
lb = [-5, -5]
ub = [5, 5]

# Equality constraint
h = (2/3)*x[0] - x[1]
include_eq = False

if include_eq:
    nlp = {'x': x, 'f': f, 'g': ca.vertcat(g1, g2, h)}
else:   
    nlp = {'x': x, 'f': f, 'g': ca.vertcat(g1, g2)}

solver_opts = {'ipopt.tol': 1e-6, 'ipopt.print_level': 0}
solver = ca.nlpsol('solver', 'ipopt', nlp, solver_opts)
x1 =  [-2, 0.8]
x2 = [-5, 0]
x3 =  [1, 0]
x4 = [1, 5]
f_opt = np.zeros(4)
x_opt = np.zeros((4, 2))

starting_points = [x1, x2, x3, x4]
for i,x0 in enumerate(starting_points):
    sol = solver(x0=x0, lbx=lb, ubx=ub)
    x_opt_current = sol['x']   
    print("Optimal solution:")
    print("x1:", x_opt_current[0])
    print("x2:", x_opt_current[1])
    print("f(x):", sol['f'])

    # Saving optimal values
    f_opt[i] = sol['f']
    x_opt[i] = x_opt_current.T


    
