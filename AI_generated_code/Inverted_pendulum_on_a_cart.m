% inverted_pendulum_lagrangian.m
% Symbolic derivation for inverted pendulum on a cart
% Coordinates: x (cart), theta (pendulum from upright, theta=0 upright)
% Parameters: M, m, l, g, F (external horizontal force on cart)

clearvars; close all; clc;
syms x theta dx dtheta ddx ddtheta real
syms M m l g F real

% Position of pendulum mass (origin at cart mass center)
xb = x + l*sin(theta);
yb = -l*cos(theta);

% Velocities
dxb = diff(xb, x)*dx + diff(xb, theta)*dtheta;   % or use jacobian
dyb = diff(yb, x)*dx + diff(yb, theta)*dtheta;

% Kinetic energy
T_cart = 1/2 * M * dx^2;
T_pend = 1/2 * m * (dxb^2 + dyb^2);
T = simplify(T_cart + T_pend);

% Potential energy: choose V(theta=0)=0 (upright zero)
V = m * g * l * (1 - cos(theta));

% Lagrangian
L = simplify(T - V);

% Generalized coordinates and velocities
q   = [x; theta];
dq  = [dx; dtheta];
ddq = [ddx; ddtheta];

% Partial derivatives for Euler-Lagrange
dLdq   = jacobian(L, q).';
dLdqdot = jacobian(L, dq).';

% Time derivatives: replace x->x(t) heuristically by using chain rule:
% d/dt(∂L/∂qdot) = ∑ ∂(∂L/∂qdot)/∂q_i * qdot_i + ∂(∂L/∂qdot)/∂qdot_i * qddot_i
% Build using Jacobians:
Dd_dt_dLdqdot = jacobian(dLdqdot, q) * dq + jacobian(dLdqdot, dq) * ddq;
EL = simplify(Dd_dt_dLdqdot - dLdq);  % equations of motion (generalized forces = [F; 0])

% Include generalized forces vector Q = [F; 0] (force applied to x only)
Q = [F; 0];
EOM = simplify(EL - Q); % should be zero

% Display simplified equations (symbolic)
eq1 = simplify(EOM(1)); % equation for x
eq2 = simplify(EOM(2)); % equation for theta

% Solve for accelerations ddx and ddtheta
sol = solve([eq1; eq2], [ddx; ddtheta], 'ReturnConditions', true);
ddx_sol = simplify(sol.ddx);
ddtheta_sol = simplify(sol.ddtheta);

% Put equations in mass matrix form: D(q)*ddq + C(q,dq) + G(q) = B*F
% Compute mass matrix D as coefficients of ddx, ddtheta in EL (before subtracting Q)
% Build symbolic matrix by collecting terms
D11 = simplify(coef(EL(1), ddx));
D12 = simplify(coef(EL(1), ddtheta));
D21 = simplify(coef(EL(2), ddx));
D22 = simplify(coef(EL(2), ddtheta));
D = [D11, D12; D21, D22];

% Compute remaining terms (move ddq terms to left: Cg = EL - D*ddq)
Cg = simplify(EL - D*ddq); % contains coriolis/centrifugal and gravity and -Q
% Separate gravity terms (depend only on q, not dq)
Gvec = simplify(subs(Cg, dtheta, 0)); % set velocities zero to isolate gravity and constant terms
% Now get velocity-dependent terms: Cvec = Cg - Gvec
Cvec = simplify(Cg - Gvec);

% Solve for mass-inertia inverse form
% Compute accelerations from D\ ( -Cvec )
acc = simplify(D \ (-Cvec + Q)); % ddq = D^{-1}*( -C + Q )
ddx_expr = simplify(acc(1));
ddtheta_expr = simplify(acc(2));

% Linearize around upright equilibrium: theta = 0, dx=0, dtheta=0
% State vector: [x; dx; theta; dtheta] -> compute A,B linearization
% Define state symbols
syms x0 dx0 theta0 dtheta0 real
% For linearization substitute symbolic variables:
subsList = {theta, dtheta, dx, x}; % note ordering
% Build state vector and dynamics: xdot = [dx; ddx_expr; dtheta; ddtheta_expr]
xdot_sym = [dx; ddx_expr; dtheta; ddtheta_expr];

% Compute Jacobians (A = df/dstate, B = df/dF)
state = [x; dx; theta; dtheta];
A = simplify(jacobian(xdot_sym, state));
B = simplify(jacobian(xdot_sym, F));

% Evaluate at equilibrium (theta=0, dtheta=0). x and dx arbitrary (don't affect linearization), set to 0.
A_lin = simplify(subs(A, {theta, dtheta, dx, x}, {0, 0, 0, 0}));
B_lin = simplify(subs(B, {theta, dtheta, dx, x}, {0, 0, 0, 0}));

% Simplify A_lin and B_lin
A_lin = simplify(A_lin);
B_lin = simplify(B_lin);

% Outputs
disp('Kinetic energy T:');
disp(T);
disp('Potential energy V:');
disp(V);
disp('Lagrangian L:');
disp(L);

disp('Equations of motion (symbolic, should equal zero):');
disp(EL);

disp('Mass matrix D(q):');
disp(simplify(D));
disp('Coriolis/Centrifugal + Other velocity terms C(q,dq):');
disp(simplify(Cvec));
disp('Gravity vector G(q):');
disp(simplify(Gvec));

disp('Accelerations (nonlinear): ddx = ');
disp(ddx_expr);
disp('ddtheta = ');
disp(ddtheta_expr);

disp('Linearized A matrix at theta=0:');
disp(A_lin);
disp('Linearized B matrix at theta=0:');
disp(B_lin);

% Helper function: coefficient extraction for symbolic expressions
function c = coef(expr, symVar)
    % Extract coefficient of symVar in expr (works for linear in symVar)
    expr = expand(expr);
    c = simplify((expr - subs(expr, symVar, 0)) / symVar);
end
