import sympy as sp
import numpy as np

# Define symbols
x, y = sp.symbols('x y')

# Define the Himmelblau function
f = (x**2 + y - 11)**2 + (x + y**2 - 7)**2

# Define the domain constraints
domain = sp.And(-5 <= x, x <= 5, -5 <= y, y <= 5)

# Compute gradient
grad_f = [sp.diff(f, var) for var in (x, y)]
x_vals = np.linspace(-5, 5, 100)
y_vals = np.linspace(-5, 5, 100)
points = [(x_val, y_val) for x_val in x_vals for y_val in y_vals if domain.subs({x: x_val, y: y_val})]

# Using nsolve to find stationary points

stationary_x, stationary_y = list(), list()

# Points represent different starting points for nsolve

for point in points:
    try:
        sol = sp.nsolve(grad_f, (x, y), point)
        curr_x, curr_y = round(sol[0].evalf(), 3), round(sol[1].evalf(), 3)

        if curr_x not in stationary_x and curr_y not in stationary_y:
            stationary_x.append(curr_x)
            stationary_y.append(curr_y)
    except sp.polys.polyerrors.PolynomialError:
        pass

stationary_points = [(x, y) for x, y in zip(stationary_x, stationary_y)]

print('Stationary points:', stationary_points)
print(len((stationary_points)))

# Classify stationary points

H = sp.hessian(f, (x, y))

classification = dict()
eigenvalues = list()

for point in stationary_points:
    H_eval = H.subs({x: point[0], y: point[1]})
    eig = H_eval.eigenvals()
    eigenvalues.append([round(i, 3) for i in eig.keys()])
    if all(val > 0 for val in eig.keys()):
        classification[point] = 'Minimum'
    elif all(val < 0 for val in eig.keys()):
        classification[point] = 'Maximum'
    else:
        classification[point] = 'Saddle'


# file for writing output to csv
output_file = 'problem4/himmelblau_extrema.csv'

classification_val = {'Maximum': 1, 'Minimum': 2, 'Saddle': 0}

with open(output_file, 'w') as file:
    file.write('x,y,Eigenvalue1,Eigenvalue2,Classification\n')
    for i, point in enumerate(stationary_points):
        c = classification[point]
        file.write(f"{point[0]},{point[1]},{eigenvalues[i][0]},{eigenvalues[i][1]},{classification_val[c]}\n")
    
print('Classification:')
for i, (point, cls) in enumerate(classification.items()):
    print(f"Point: {point}, Eigenvalues: {eigenvalues[i]}, Class: {cls}\n")





