
# 1 - Lagrange for General Box-Constrained NLP
General NLP
$$ \begin{align*}
&\min_{x} f(x) \\\\
&\text{s.t.} \quad h(x) = 0,\\
&\quad g_{l} \leq g(x) \leq g_{u},\\
&\quad x_{l} \leq x \leq x_{u}\\
 \end{align*}$$
Definition of Lagrangian
$$\mathcal{L}(x, \lambda) = f(x) - \sum_{i \in \mathcal{I} \cup \mathcal{E}} \lambda_i c_i(x)$$
Where $\forall i \in \mathcal{I}:$ $c_{i}(x) \geq 0$.
## Dimensions
- $x \in \mathbb{R}^n$
- $f: \mathbb{R}^{n}\mapsto \mathbb{R} \quad f \in \mathcal{C}^2$
- $h: \mathbb{R}^{n}\mapsto \mathbb{R}^{m_{\mathcal{E}}} \quad h \in \mathcal{C}^2$
-  $g: \mathbb{R}^{n}\mapsto \mathbb{R}^{m_{\mathcal{I}}} \quad h \in \mathcal{C}^2$
($\cdot \in \mathcal{C}^2$ as a result of smoothness necessary for SQP)
We therefore convert the nonlinear inequality and box constraints into this form:
### Equality Constraints

$$c_{i \in \mathcal{E}}(x) = h(x) = 0$$
### Inequality Constraint

$$\begin{align*}
g_{l} \leq g(x) &\implies g(x) - g_{l} \geq 0 \\
g_{u} \geq g(x) &\implies  g_{u} - g(x) \geq 0
\end{align*}$$
$$
 \begin{align*} g(x) - g_{l} &\geq 0 \\
g_{u} - g(x) &\geq 0\\
 \end{align*}
$$
Expressing as a vector
$$
\begin{align*}
c_{i \in \mathcal{I}}(x) &\geq 0 \\
\begin{bmatrix} g(x) - g_{l} \\ g_{u} - g(x) \end{bmatrix} &\geq 0
\end{align*}
$$
Where
$$c_{i \in \mathcal{I}}(x) =
\begin{bmatrix} g(x) - g_{l} \\ g_{u} - g(x) \end{bmatrix} \in \mathbb{R}^{m_{\mathcal{I}}} $$
### Box Constraints
$$
\begin{align*}
x - x_{l} &\geq 0\\
x_{u} - x &\geq 0
\end{align*}
$$
## Lagrange
### Multipliers

- $\mathrm{z} \in \mathbb{R}^{m_{\mathcal{I}}}$: Inequality constraints.
- $\mathrm{y} \in \mathbb{R}^{m_\mathcal{E}}$: Equality constraints.
- $\mathrm{a} \in \mathbb{R}^{n}$: Lower box constraints.
- $\mathrm{b} \in \mathbb{R}^{n}$: Upper box constraint.
$$ 
\begin{align*}
\mathcal{L}(x, \mathrm{y}, \mathrm{z}, \mathrm{a}, \mathrm{b}) &= f(x) - \mathrm{y}'c_{\mathcal{E}}(x)-  \mathrm{z}' c_{\mathcal{I}}(x) - \mathrm{a}'(x-x_{l}) -\mathrm{b}'(x_{u} -x)
\end{align*}
$$
Where $c_{\mathcal{E}}(x) = h(x)$ and $c_{i \in \mathcal{I}}(x) =\begin{bmatrix} g(x) - g_{l} \\ g_{u} - g(x) \end{bmatrix}$. 
# 2 - Necessary First-Order Optimality Conditions
# 3 - Sufficient Second-Order Optimality Conditions
# 4 - Himmelblau's test problem 

We use the example outlined in Lecture 2C, along with a linear equality constraint, passing through one of the minima:

$$\begin{align*}
&\min_{x_{1}, x_{2}} f(x_{1}, x_{2}) = (x_{1}^{2}+x_{2}-11)^{2}+(x_{1}+x_{2}^{2}-7)^{2} 
\end{align*}$$
s.t.
$$
\begin{align*}
h(x)=\frac{2}{3}x_{1}-x_{2}&=0\\
c_{1}(x, y) = (x_{1}+2)^{2}-x_{2} &\geq 0 \\
c_{2}(x, y) = -4x_{1}+10x_{2} &\geq 0 \\
\begin{bmatrix}-5 \\ -5 \end{bmatrix} \leq \begin{bmatrix}x_{1}\\x_{2}\end{bmatrix} \leq \begin{bmatrix}5 \\ 5\end{bmatrix}
\end{align*}
$$

## Re-Writing
Writing this in the standard form defined above, we have
- $h = \frac{2}{3}x_{1}-x_{2}$
- $g(x) = \begin{bmatrix} (x_{1}+2)^{2}-x_{2} \\ -4x_{1}+10 x_{2}\end{bmatrix}$
- $x_{l} = \begin{bmatrix}-5 & -5\end{bmatrix}'$
- $x_{u}= \begin{bmatrix}5 & 5\end{bmatrix}'$
- $g_{l}= \begin{bmatrix}0 & 0\end{bmatrix}'$
- $g_{u} = \begin{bmatrix} \infty & \infty\end{bmatrix}'$ (Note that in practice this is simply set to some very large value)
s.t.
$$\frac{2}{3}x_{1}-x_{2}=0$$
$$
\begin{align*}
c_{\mathcal{I}}(x) = G(x) = 
\begin{bmatrix} g(x) - g_{l} \\ g_{u} - g(x) \end{bmatrix} = &\geq 0 \rightarrow \begin{bmatrix} (x_{1}+2)^{2}-x_{2} \\
 -4x_{1}+10 x_{2} \\
 \infty - ((x_{1}+2)^{2}-x_{2}) \\
 \infty - (-4x_{1}+10 x_{2}) \end{bmatrix} \geq 0
\end{align*}
$$
## Lagrange
$$ \mathcal{L}(x, \mathrm{y}, \mathrm{z}, \mathrm{a}, \mathrm{b}) = (x_{1}^{2}+x_{2}-11)^{2}+(x_{1}+x_{2}^{2}-7)^{2} - \mathrm{y}' \cdot \begin{bmatrix} \frac{2}{3}x_{1}-x_{2}\end{bmatrix} - \mathrm{z}' \begin{bmatrix} (x_{1}+2)^{2}-x_{2} \\
 -4x_{1}+10 x_{2} \\
 \infty - ((x_{1}+2)^{2}-x_{2}) \\
 \infty - (-4x_{1}+10 x_{2}) \end{bmatrix}  - \mathrm{a}'\begin{bmatrix}x_{1} + 5 \\ x_{2} + 5\end{bmatrix} - \mathrm{b}' \begin{bmatrix}5 - x_{1} \\ 5 - x_{2}\end{bmatrix} $$
### Gradient
$$
\begin{align*}
\nabla \mathcal{L}(x, \mathrm{y}, \mathrm{z}, \mathrm{a}, \mathrm{b}) = \begin{bmatrix} 4 x_{1} (x_{1}^{2}+x_{2} -11) \\ 2 (x_{1}^{2}+x_{2} -11)\end{bmatrix} + \begin{bmatrix} 2(x_{1}+x_{2}^{2}-7) \\ 4 x_{2} (x_{1}+x_{2}^{2}-7\end{bmatrix} - \begin{bmatrix}\frac{2}{3} \\ -1\end{bmatrix} \mathrm{y}- \begin{bmatrix}2 x_{1} +4 & -1 \\ -4 & 10 \\ -2 x_{1}-4 & 1 \\ 4 & 10\end{bmatrix}^{\top} \mathrm{z} - \mathrm{a} + \mathrm{b}
\end{align*}
$$
### Hessian (of the Lagrange)
Using the definition
$$\nabla_{x x}^{2} \mathcal{L}(x, \mathrm{y}, \mathrm{z}, \mathrm{a}, \mathrm{b}) = \nabla^{2} f(x) - \sum\limits_{i=1}^{m_\mathcal{E}} \mathrm{y}_{i} \nabla^{2}h_{i}(x) - \sum\limits_{i=1}^{m_{\mathcal{I}}} \mathrm{z}_{i} \nabla^{2}G_{i}(x)$$
In our case:
- $m_{\mathcal{I}} = 4$
- $m_{\mathcal{E}}=1$
- $\nabla^{2}h(x) = \begin{bmatrix}0 \\ 0\end{bmatrix}$
- $\forall i \in m_\mathcal{I}: \ \nabla^{2}G_{i}(x)$ are $2 \times 2$ "block matrices"
	- $\nabla^{2}G_{1}(x) = \begin{bmatrix}2 & 0 \\ 0 & 0\end{bmatrix}$
	- $\nabla^{2}G_{2}(x) = 0$ (all entries 0)
	- $\nabla^{2}G_{3}(x) = \begin{bmatrix}-2 & 0 \\ 0 & 0\end{bmatrix}$
	- $\nabla^{2}G_{4}(x)=0$

## Contour with Stationary points
(Done)
# 5 - Solution using fmincon
(Done)
- Check Lagrangian hessian (42:30 in lecture 09A) - `fmincon` has sign of Lagrange multipliers swapped.
## Options

## Algorithms
Saving the results of all algorithms and their iterates
### SQP 
- No Hessian can be passed - Must make BFGS approximation.
### Interior Point (default)
### Active Set
### Trust-region reflective
# 6 - SQP with line search
## Theory
- SQP is a Newton's method on the KKT-conditions for the nonlinear program:
$$x_{n+1} = x_{n} - \frac{f(x_{n})}{f'(x_{n})}$$
- Each step is a QP (rather than a linear system of equations/unconstrained QP)
- Major step in SQP by the quasi-Newton method is a solution of convex QP with equality and inequality constraints. 
	- Need to enforce convexity of QP.
	- Solution of this QP involves sequence of equality constrained convex QPs (active set etc.)
- Relies on linearised system.
	- Line search helps with convergence.
	- Line search **guarantees** convergence.
"Big picture"
$$\text{Nonlinear problem} \rightarrow \text{Inequality Constrained QPs} \rightarrow \text{Equality Constrained QPs} \rightarrow \text{KKT solutions}$$
### Problem set up 
$$ \begin{align*}
&\min_{x} f(x) \\\\
&\text{s.t.} \quad h(x) = 0,\\
&\quad g_{l} \leq g(x) \leq g_{u},\\
&\quad x_{l} \leq x \leq x_{u}\\
 \end{align*}$$
- $x_{l}, x_{u}, x \in \mathbb{R}^n$
- $f: \mathbb{R}^{n}\mapsto \mathbb{R} \quad f \in \mathcal{C}^2$
- $h: \mathbb{R}^{n}\mapsto \mathbb{R}^{m_\mathcal{E}} \quad h \in \mathcal{C}^2$
-  $g: \mathbb{R}^{n}\mapsto \mathbb{R}^{m_{\text{NL box}}} \quad h \in \mathcal{C}^2$
- $g_{l}, g_{u} \in \mathbb{R}^{m_{\text{NL box}}}$

#### A Comment on the Dimensionality of $g$ and the Associated Lagrangian Multiplier $\mathrm{z}$

We define the number of inequality constraints as $m_{\text{Nonlinear Box}}$. When this is converted to regular inequality constraints, we end up with a doubling of the number of inequality constraints, compared to the number of box constraints.

$$g_{l} \leq g(x) \leq g_{u} \rightarrow \begin{bmatrix} g(x) - g_{l} \\ g_{u} - g(x) \end{bmatrix} \geq 0$$
As a consequence, the dimensions of the **Lagrange Multiplier** associated with the nonlinear inequality constraints $\mathrm{z} \in \mathbb{R}^{m_{\mathcal{I}}}$ where $m_\mathcal{I} = 2 \cdot m_{\text{NL box}}$

### Writing in Standard Form
Standard form of a box-constrained NLP (equality and inequality constraints)

$$
\begin{align*}
&\min_{x \in \mathbb{R}^{n}} f(x) \\
&\text{s.t.} \quad h(x) = 0 \\
&\quad \quad g(x) \geq 0 \\
&\quad x_{l} \leq x \leq x_{u}
\end{align*}
$$
Converting for our problem:
- $h$ is unchanged.
- $g$ is altered to form two inequality constraints.
- Box constraints form two more inequality constraints.
$$\begin{align*}
c_{i \in \mathcal{I}}(x) =
\begin{bmatrix} g(x) - g_{l} \\ g_{u} - g(x) \end{bmatrix} &\geq 0\\
x - x_{l} &\geq 0\\
x_{u} - x &\geq 0
\end{align*}$$
Final form

$$
\min_{x \in \mathbb{R}^{n}} f(x)
$$
Subject to
$$
\begin{align*}
h(x) &= 0 \\
c_{i \in \mathcal{I}}(x) =
\begin{bmatrix} g(x) - g_{l} \\ g_{u} - g(x) \end{bmatrix} &\geq 0\\
x - x_{l} &\geq 0\\
x_{u} - x &\geq 0
\end{align*}
$$
### Lagrangian

$$ 
\begin{align*}
\mathcal{L}(x, \mathrm{z}, \mathrm{z}, \mathrm{a}, \mathrm{b}) &= f(x) - \mathrm{\mathrm{z}}' c_{\mathcal{E}}(x) - \mathrm{\mathrm{z}}'c_{\mathcal{I}}(x) - \mathrm{a}'(x - x_{l}) - \mathrm{b}'(x_{u} - x)\\
&= f(x) - \mathrm{\mathrm{y}}' h(x) - \mathrm{\mathrm{z}}'\begin{bmatrix} g(x) - g_{l} \\ g_{u} - g(x) \end{bmatrix} - \mathrm{a}'(x - x_{l}) - \mathrm{b}'(x_{u} - x)
\end{align*}$$
Where the following **Lagrange multipliers** correspond to
- $\mathrm{y} \in \mathbb{R}^{m_\mathcal{E}}$: Equality constraints.
- $\mathrm{z} \in \mathbb{R}^{m_\mathcal{I}}$: Inequality constraints.
- $\mathrm{a} \in \mathbb{R}^n$: Lower bound.
- $\mathrm{b} \in \mathbb{R}^n$: Upper bound.
### KKT conditions (first-order necessary)
**Stationarity:**
$$ 
\begin{align*}
\nabla_{x} \mathcal{L}(x, \mathrm{y}, \mathrm{z}, \mathrm{a}, \mathrm{b}) &= \nabla f(x) -  \nabla h(x) \mathrm{\mathrm{y}} -  \nabla\begin{bmatrix} g(x) - g_{l} \\ g_{u} - g(x) \end{bmatrix} \mathrm{\mathrm{z}} -  \nabla(x - x_{l}) \mathrm{a}-  \nabla(x_{u} - x) \mathrm{b} \\
&= \nabla f(x) - \nabla h(x)\mathrm{y} - \begin{bmatrix}  \nabla g(x) \\- \nabla g(x) \end{bmatrix} \mathrm{z} - \mathrm{a} + \mathrm{b}
\end{align*}$$
Checking dimensions:
- $\nabla f(x) \in \mathbb{R}^{n \times 1}$
-  $\mathrm{y} \in \mathbb{R}^{m_{\mathcal{E}} \times 1}$, $\nabla h(x) \in \mathbb{R}^{m_{\mathcal{E}} \times n} \implies \nabla h(x) \mathrm{y} \in \mathbb{R}^{n \times 1}$  
- $\mathrm{z} \in \mathbb{R}^{m_{\mathcal{I}} \times 1}$, $\begin{bmatrix} \nabla g(x)\\- \nabla g(x) \end{bmatrix} \in \mathbb{R}^{m_{\mathcal{I}} \times n} \implies  \begin{bmatrix} \nabla g(x)\\- \nabla g(x) \end{bmatrix} \mathrm{z} \in \mathbb{R}^{n \times 1}$
- $\mathrm{a}, \mathrm{b} \in \mathbb{R}^{n \times 1}$ 
**Primal Feasibility:**
$$\begin{align*}
h(x) &= 0\\
\begin{bmatrix} g(x) - g_{l} \\ g_{u} - g(x) \end{bmatrix} &\geq 0 \\
x - x_{l} &\geq 0\\
x_{u} - x &\geq 0
\end{align*}$$
**Dual Feasibility:**
$$\mathrm{z}_{i} \geq 0, \quad i = 1, \dots, m_{\mathcal{I}}$$
**Complementary Slackness:**
$$
\begin{align*}
\mathrm{z}_{i} = 0 \lor c_i(x) = 0 \iff \mathrm{z}_i c_i (x) &= 0, \quad i = 1, \dots, m_{\mathcal{I}}\\
\mathrm{a}_{i} (x_{i} - x_{l, i}) &= 0, \quad i = 1, \dots, n\\
\mathrm{b}_{i} (x_{u, i} - x_{i}) &= 0, \quad i = 1, \dots, n
\end{align*}
$$
## Practical Implementation
### The problem
$$
\min_{x \in \mathbb{R}^{n}} f(x)
$$
Subject to
$$
\begin{align*}
h(x) &= 0 \\
c_{\mathcal{I}}(x) = G(x) = 
\begin{bmatrix} g(x) - g_{l} \\ g_{u} - g(x) \end{bmatrix} &\geq 0\\
x - x_{l} &\geq 0\\
x_{u} - x &\geq 0
\end{align*}
$$
### Lagrangian

$$
\begin{align*}
\mathcal{L}(x, \mathrm{y}, \mathrm{z}, \mathrm{a}, \mathrm{b}) &= f(x) - \mathrm{\mathrm{y}}' h(x) - \mathrm{\mathrm{z}}'\begin{bmatrix} g(x) - g_{l} \\ g_{u} - g(x) \end{bmatrix} - \mathrm{a}'(x - x_{l}) - \mathrm{b}'(x_{u} - x)\\
\nabla_{x} \mathcal{L}(x, \mathrm{y}, \mathrm{z}, \mathrm{a}, \mathrm{b}) &= \nabla f(x) - \nabla h(x)\mathrm{y} - \begin{bmatrix}  \nabla g(x) \\- \nabla g(x) \end{bmatrix} \mathrm{z} - \mathrm{a} + \mathrm{b} \\
\nabla_{x x}^{2} \mathcal{L}(x, \mathrm{y}, \mathrm{z}, \mathrm{a}, \mathrm{b}) &= \nabla^{2} f(x) - \sum\limits_{i=1}^{m_\mathcal{E}} \mathrm{y}_{i} \nabla^{2}h_{i}(x) - \sum\limits_{i=1}^{m_{\mathcal{I}}} \mathrm{z}_{i} \nabla^{2}G_{i}(x)
\end{align*}
$$
### Quadratic Program - Finding Search Direction
Found via the solution of the following QP: (18.11 N&W)

$$\min_{p_k \in \mathbb{R}^{n}} \quad \frac{1}{2} p_k' \left[\nabla_{x x}^{2} \mathcal{L}(x, \mathrm{y}, \mathrm{z}, \mathrm{a}, \mathrm{b})\right] p_k + \left[\nabla f(x^{k})\right]p_k + f(x^k)$$
s.t.
$$
\begin{align*}
\left[\nabla h(x^{k})\right]' p_k + h(x^{k}) &= 0\\
\left[ \nabla G(x^{k})\right]' p_k + G(x^{k}) &\geq 0
\end{align*}
$$
Allows us to find the optimal step direction in $x$ ($p_k$) and updated Lagrange multipliers.

#### $\ell^{1}$ Quadratic Sub-Problem
In order to enforce feasibility, we reformulate the quadratic subproblem as a an $\ell^{1}$ penalty problem
$$\min_{p_{k},v,w,t} \frac{1}{2} p_{k}'\nabla_{xx}^{2}\mathcal{L}p_{k} + \nabla f_{k}'p_{k} + \mu(e_{v}' v +e_{w}'w+e_{t}'t)$$
s.t.
$$\begin{align*}
\left[\nabla h(x_{k}\right]'p_{k}+h(x_{k})&=v-w\\
\left[\nabla G(x_{k}\right]'p_{k} + G(x_{k}) &\geq -t\\
v,w,t &\geq 0
\end{align*}$$

### Quasi-Newton Method - BFGS and Modified BFGS Update for Hessian
Will be used in second method for approximating Hessian.
$$
\begin{align*}
p &= x^{k+1} - x^{k} = \alpha \cdot p_{k}\\
q &= \nabla_{x} \mathcal{L}(x^{k+1}, \mathrm{y}^{k+1}, \mathrm{z}^{k+1}, \mathrm{a}^{k+1}, \mathrm{b}^{k+1}) - \nabla_{x} \mathcal{L}(x^{k}, \mathrm{y}^{k+1}, \mathrm{z}^{k+1}, \mathrm{a}^{k+1}, \mathrm{b}^{k+1})
\end{align*}
$$
#### BFGS 
$$B \leftarrow B + \frac{qq'}{p'q}-\frac{(Bp)(Bp)'}{p'(Bp)}$$
#### Modified BFGS Update
$$\begin{align*}
\theta &= \begin{cases} 1 \quad \text{if} \quad p'q \geq 0.2 p'(Bp)\\
\dfrac{0.8p'(Bp)}{p'(Bp)-p'q} \quad \text{otherwise}\end{cases}\\
r &= \theta q + (1-\theta)(Bp) \\
B &\leftarrow B + \frac{rr'}{p'r} - \frac{(Bp)(Bp)'}{p'(Bp)}
\end{align*}$$
### Local SQP - Pseudocode

### Line Search - Computing Optimal Step Direction

New iterate (let *step length* be $a$)

$$x = x^{k}+a p_k^{k}$$
Merit Function (**Powell's $\ell_{1}$**)
$$\begin{align*}
\phi(\alpha) = P(x, \mathrm{y}, \mathrm{z}, \mathrm{a}, \mathrm{b}) &= P(x^{k}+ \alpha p_k^{k}, \mathrm{y}, \mathrm{z}, \mathrm{a}, \mathrm{b})\\
&= f(x^{k}+ \alpha p_k^{k})+\lambda' \vert h(x^{k}+\alpha p_k^{k})\vert + \mu' \vert \min{\{0, c_{\mathcal{I}}(\mathrm{x^{k}+\alpha p_k^{k}})\}} \vert
\end{align*} $$
With
- $\lambda \geq \vert \mathrm{y} \vert$
- $\mu \geq \vert z \vert$
**Powell's update** of penalty parameters
$$\begin{align*}
\lambda &= \max \left\{\vert \mathrm{y} \vert, \frac{1}{2}(\lambda + \vert \mathrm{y} \vert) \right\}\\
\mu &= \max \left\{\vert z\vert, \frac{1}{2}(\mu+\vert \mathrm{z}\vert)\right\}
\end{align*}$$

### Overview
1. **Initialization**:
   - Choose an initial guess $x_{0}$ for the optimisation variables.
   - Initialise Lagrange multipliers, penalty parameters, and other algorithm parameters.
2. **Main Iteration Loop**:
   - **Step 1: Compute Search Direction**:
     - Linearise the constraints and the objective function around the current iterate $x_{k}$.
     - Formulate and solve a Quadratic Programming (QP) subproblem to obtain the search direction $p_k$
   - **Step 2: Line Search**:
     - Perform a line search along the search direction $p_{k}$ to determine the step size $\alpha_{k}$.
     - Update the iterate:  $x_{k+1} = x_k + \alpha_k p_k$ .
   - **Step 3: Update Lagrange Multipliers**:
     - Update Lagrange multipliers associated with the active constraints based on the new iterate $x_{k+1}$ .
   - **Step 4: Termination Check**:
     - Check for convergence criteria, such as small changes in the objective function, feasibility of constraints, or reaching a maximum number of iterations.
     - If convergence criteria are met, terminate the algorithm; otherwise, return to Step 1.
3. **Necessary Sub-algorithms**:
   - **Quadratic Programming (QP) Solver**:
     - Solves the QP subproblem in each iteration to obtain the search direction.
   - **Line Search Algorithm**:
     - Determines the step size along the search direction to ensure sufficient decrease in the objective function.
   - **Constraint Handling**:
     - Ensures that constraints are satisfied and handles active and inactive constraints appropriately.
   - **Convergence Criteria**:
     - Defines conditions for terminating the algorithm, such as reaching a specified tolerance level for optimality and feasibility, or achieving a maximum number of iterations.
4. **Termination**:
   - The algorithm terminates when convergence criteria are satisfied.

The SQP algorithm iteratively improves the solution by approximating the original nonlinear optimisation problem with a sequence of QP subproblems. Each subproblem accounts for the nonlinear constraints and the objective function, and the solution gradually converges to a locally optimal solution. Proper handling of constraints, choice of QP solver, line search strategy, and convergence criteria are essential for the effectiveness and efficiency of the algorithm.
### Parameters
- objective function: $f(x)$.
- Initial guess: $x_0$.
- lower and upper bounds (box constraint).
- nonlinear upper and lower bounds.
- Equality constraints.
- Hessian (user defined hessian or BFGS) (optional - default is BFGS)
- Tolerance (optional)
- Maximum number of iterations (optional)
### Output
- `x`: optimal solution
- `fval`: optimal function value
- `exitflag`: exit flag
	- 1: converged
	- 0: maximum number of iterations reached
- `iterations`: number of iterations
- `lambda`: Struct containing Lagrange multipliers
	- `lambda.ineq`: Lagrange multipliers for inequality constraints
	- `lambda.eq`: Lagrange multipliers for equality constraints
	- `lambda.lower`: Lagrange multipliers for lower bounds
	- `lambda.upper`: Lagrange multipliers for upper bounds
- `grad`: gradient at the solution
- `hessian`: hessian at the solution
## Pseudocode
Handle inputs/optional parameters.
- $\text{tol} = 10^{-6}$
- $\text{MaxIter} = 100$
- Default Hessian is a (damped) *BFGS* approximation 
- If given $x_0$ is not feasible, place in centre of box constraints.
#### Initialise Variables
##### Information
- $\text{exitflag = 0}$ (has not converged)
- $\text{iterations} = 0$ 
- $\text{stop} = \text{False}$
##### Empty Arrays
- $x, \nabla f \in \mathbb{R}^{n \times \text{MaxIter}}$: Storing $x$ and gradient values.
- $f \in \mathbb{R}^{1 \times \text{MaxIter}}$: Storing objective function values.
##### Initial Values
- $\lambda =0$ for all Lagrange multipliers (Equality, inequality, lower, upper)
- $h_{0}, \nabla h_{0} = h(x_{0}), \nabla h(x_{0})$ 
- $G_{0}, \nabla G_{0} = G(x_{0}), \nabla G(x_{0})$
#### Solving Local QP
##### With quadprog
Form of QP for quadprog:
$$\min_{x \in \mathbb{R}^{n}} \frac{1}{2}x' H x + f' x$$
s.t.
$$\begin{align*}
A \cdot x &\leq \mathrm{b} \\
Aeq \cdot x &= \mathrm{beq}\\
lb \leq x &\leq ub
\end{align*}$$
For SQP:

$$\min_{p_k \in \mathbb{R}^{n}} \quad \frac{1}{2} p_k' \left[\nabla_{x x}^{2} \mathcal{L}(x, \mathrm{y}, \mathrm{z})\right] p_k + \left[\nabla f(x^{k})\right]p_k + f(x^k)$$
s.t:  **wrong!**
$$\begin{align*}
-\left( \left[ \nabla G(x^{k})\right]' p_k + G(x^{k})\right) \leq 0 &\implies  -\left[ \nabla G(x^{k})\right]' p_k ) \geq G(x^{k})\\
\left[\nabla h(x^{k})\right]' p_k + h(x^{k}) = 0 &\implies \left[\nabla h(x^{k})\right]' p_{k} = -h(x^{k}) \\
x_{l} \leq x \leq x_{u} &\implies x - x_{l} \leq p_{k} \leq x_{u}- x
\end{align*}$$
quadprog inputs and their associated variables:
**Inequality constraint:**
- $A \rightarrow - \nabla G(x^{k})'$
- $\mathrm{b} \rightarrow G(x^{k})$
**Equality Constraints:**
- $Aeq \rightarrow \left[\nabla h(x^{k})\right]'$
- $\mathrm{beq} \rightarrow - h(x^{k})$
### Update step
#### Line-Search


# 7 - Trust Region SQP
