clear all;

syms F
% p: position, v: velocity, theta: angle, omega: angular velocity
% F: input force

M = 0.5;
m = 0.2;
l = 0.3;
I = 0.006;
g = 9.81;

q = [sym('x'); sym('theta')];
q_dot = [sym('x_dot'); sym('theta_dot')];
q_ddot = [sym('x_ddot'); sym('theta_ddot')];

% 1. Kinematics
% Position of the pendulum mass (x, y)
x_p = q(1) + l*sin(q(2));
y_p = l*cos(q(2));

% Velocities (time derivatives)
vx_p = diff(x_p, q(1))*q_dot(1) + diff(x_p, q(2))*q_dot(2);
vy_p = diff(y_p, q(1))*q_dot(1) + diff(y_p, q(2))*q_dot(2);

% 2. Energy Equations
T_cart = 0.5 * M * q_dot(1)^2;
T_pend_trans = 0.5 * m * (vx_p^2 + vy_p^2);
T_pend_rot = 0.5 * I * q_dot(2)^2;
T = T_cart + T_pend_trans + T_pend_rot;

V = m * g * y_p;

L = T - V;

% 3. Euler-Lagrange Equations
% State vector q = [p; theta], q_dot = [v; omega]

% Partial derivatives for Lagrange
dL_dqdot = jacobian(L, q_dot)';
dL_dq = jacobian(L, q)';

% Time derivative of (dL/dq_dot)
% Since dL_dqdot depends on p, v, theta, and omega:
dt_dL_dqdot = jacobian(dL_dqdot, [q; q_dot]) * [q_dot; q_ddot];

% Equations of motion: d/dt(dL/dq_dot) - dL/dq = Q_external
% Q_external is [F; 0] (Force applied to cart)
EOM = dt_dL_dqdot - dL_dq - [F; 0];

% Solve for accelerations [v_dot; omega_dot]
sol = solve(EOM == 0, q_ddot);

% Define state-space: x_dot = f(x, u)
f = [q_dot(1); sol.x_ddot; q_dot(2); sol.theta_ddot];

% Compute Jacobians for A and B
states = [q(1); q_dot(1); q(2); q_dot(2)];
A_sym = jacobian(f, states)
B_sym = jacobian(f, F);

% Substitute physical values
params = {q(1), q_dot(1), q(2), q_dot(2)};
values = {0, 0, 0, 0}; % Linearize at theta=0

A = double(subs(A_sym, params, values));
B = double(subs(B_sym, params, values));

% LQR Design
Q = diag([10, 1, 100, 1]); % Penalize position and angle
R = 0.1;
K = lqr(A, B, Q, R);

disp('LQR Gain Matrix K:');
disp(K);