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

% --- Inputs assumed defined earlier ---
% q   : n×1 vector of generalized coordinates symbolic
% dq  : n×1 vector of generalized velocities symbolic
% ddq : n×1 vector of generalized accelerations symbolic
% EL  : n×1 Euler-Lagrange vector (should equal Q, i.e. EL - Q = 0)
% Q   : n×1 generalized forces vector (e.g., [F; 0; ...])
% Example for pendulum: q = [x; theta]; dq = [dx; dtheta]; ddq = [ddx; ddtheta];

n = numel(q);

% Mass matrix D(q): coefficients of accelerations in EL (linear in ddq)
D = simplify(jacobian(EL, ddq));  % n x n

% Remaining terms: move D*ddq to left, Cg contains velocity and gravity and -Q
Cg = simplify(EL - D*ddq);        % n x 1

% Gravity vector G(q): set velocities to zero to isolate q-only terms
zeroDQ = sym(zeros(size(dq)));
Gvec = simplify(subs(Cg, dq, zeroDQ));  % n x 1

% Velocity-dependent terms C(q,dq) (Coriolis/centrifugal + other velocity terms)
Cvec = simplify(Cg - Gvec);  % n x 1

% Solve for nonlinear accelerations (ddq = D^{-1}*( -C + Q ))
% Ensure Q is available and consistent
acc_nl = simplify(D \ (-Cvec + Q));   % n x 1 symbolic ddq expressions

% Linearization about equilibrium (q0, dq0). Use symbolic q0,dq0 or numeric later.
% Define equilibrium point here (example upright at zero)
q0  = sym(zeros(size(q)));    % change if needed
dq0 = sym(zeros(size(dq)));

% Evaluate D at equilibrium
D0 = simplify(subs(D, q, q0));  % D evaluated at q0 (no dq dependence)

% Linearize C: keep first-order in dq -> C_lin = (∂C/∂dq)|0 * dq
C_dq_jac = simplify(jacobian(Cvec, dq));   % n x n
C_lin = simplify(C_dq_jac * dq);           % n x 1 (first-order in dq)

% Linearize G: G_lin = (∂G/∂q)|0 * (q - q0)
G_q_jac = simplify(jacobian(Gvec, q));    % n x n
% use small displacement delta_q = q - q0; here we keep symbolic q and assume q0 zeros
delta_q = q - q0;
G_lin = simplify(G_q_jac * delta_q);      % n x 1 (first-order in q)

% Solve for mass-inertia inverse form
% Compute accelerations from D\ ( -Cvec )
acc = simplify(D \ (-Cvec + Q)); % ddq = D^{-1}*( -C + Q )
ddx_expr = simplify(acc(1));
ddtheta_expr = simplify(acc(2));

% state: q = [x; theta], dq = [dx; dtheta]
Bvec = [1; 0];

% Linearized implicit dynamics: D0*ddq + C_lin + G_lin = Q_lin
% For input linearization, linearize Q if needed (here assume Q linear in input F: Q = Bvec*F)
% If Q depends on q or dq include jacobian terms similarly. For simplicity assume Q = Bvec*F:
% Bvec must be n x m (m = number of inputs). Example: Bvec = [1;0] for single force.
% Ensure Bvec is defined
% Solve for linear ddq: ddq = D0 \ ( -C_lin - G_lin + Bvec*F )
ddq_lin = simplify(D0 \ ( -C_lin - G_lin + Bvec*F ));  % n x 1 (affine in q, dq, F)

% Build state vector and linear state-space (for mechanical systems typical state = [q; dq])
state = [q; dq];              % 2n x 1
% xdot linear: [dq; ddq_lin]
xdot_lin = [dq; ddq_lin];    % 2n x 1

% Compute A,B matrices symbolically
A_lin = simplify(jacobian(xdot_lin, state));   % 2n x 2n
B_lin = simplify(jacobian(xdot_lin, symvar(Bvec*F))); % 2n x m, more robust below

% Simpler: jacobian wrt inputs (assumes inputs are named symbols appearing in Bvec*F)
B_lin = simplify(jacobian(xdot_lin, F));       % 2n x m (F may be scalar or vector)

% Evaluate A,B at equilibrium (substitute q->q0, dq->dq0)
A_lin = simplify(subs(A_lin, [q; dq], [q0; dq0]));
B_lin = simplify(subs(B_lin, [q; dq], [q0; dq0]));

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

disp('Nonlinear accelerations ddq = '); disp(acc_nl);

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
