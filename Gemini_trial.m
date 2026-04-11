syms p v theta omega F M m l I g real
% p: position, v: velocity, theta: angle, omega: angular velocity
% F: input force

% 1. Kinematics
% Position of the pendulum mass (x, y)
x_p = p + l*sin(theta);
y_p = l*cos(theta);

% Velocities (time derivatives)
vx_p = diff(x_p, p)*v + diff(x_p, theta)*omega;
vy_p = diff(y_p, p)*v + diff(y_p, theta)*omega;

% 2. Energy Equations
T_cart = 0.5 * M * v^2;
T_pend_trans = 0.5 * m * (vx_p^2 + vy_p^2);
T_pend_rot = 0.5 * I * omega^2;
T = T_cart + T_pend_trans + T_pend_rot;

V = m * g * y_p;

L = T - V;

% 3. Euler-Lagrange Equations
% State vector q = [p; theta], q_dot = [v; omega]
q = [p; theta];
q_dot = [v; omega];

% Partial derivatives for Lagrange
dL_dqdot = [diff(L, v); diff(L, omega)];
dL_dq = [diff(L, p); diff(L, theta)];

% Time derivative of (dL/dq_dot)
% Since dL_dqdot depends on p, v, theta, and omega:
dt_dL_dqdot = jacobian(dL_dqdot, [q; q_dot]) * [v; omega; sym('v_dot'); sym('omega_dot')];

% Equations of motion: d/dt(dL/dq_dot) - dL/dq = Q_external
% Q_external is [F; 0] (Force applied to cart)
EOM = dt_dL_dqdot - dL_dq - [F; 0];

% Solve for accelerations [v_dot; omega_dot]
vars = [sym('v_dot'); sym('omega_dot')];
sol = solve(EOM == 0, vars);

% Define state-space: x_dot = f(x, u)
f = [v; sol.v_dot; omega; sol.omega_dot];

% Compute Jacobians for A and B
states = [p; v; theta; omega];
A_sym = jacobian(f, states);
B_sym = jacobian(f, F);

% Substitute physical values
params = {M, m, l, I, g, p, v, theta, omega};
values = {0.5, 0.2, 0.3, 0.006, 9.81, 0, 0, 0, 0}; % Linearize at theta=0

A = double(subs(A_sym, params, values));
B = double(subs(B_sym, params, values));

% LQR Design
Q = diag([10, 1, 100, 1]); % Penalize position and angle
R = 0.1;
K = lqr(A, B, Q, R);

disp('LQR Gain Matrix K:');
disp(K);