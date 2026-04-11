syms p v theta omega F1 F2 M m l I g real
% F1: Force on cart
% F2: Torque at pendulum joint

% --- 1. Kinematics & Energy (Same as before) ---
x_p = p + l*sin(theta);
y_p = l*cos(theta);
vx_p = diff(x_p, p)*v + diff(x_p, theta)*omega;
vy_p = diff(y_p, p)*v + diff(y_p, theta)*omega;

T = 0.5*M*v^2 + 0.5*m*(vx_p^2 + vy_p^2) + 0.5*I*omega^2;
V = m*g*y_p;
L = T - V;

% --- 2. Euler-Lagrange with Two Inputs ---
q = [p; theta];
q_dot = [v; omega];
dL_dqdot = [diff(L, v); diff(L, omega)];
dL_dq = [diff(L, p); diff(L, theta)];

% Time derivative
dt_dL_dqdot = jacobian(dL_dqdot, [q; q_dot]) * [v; omega; sym('v_dot'); sym('omega_dot')];

% External Forces: F1 acts on p, F2 acts on theta
Q_ext = [F1; F2]; 
EOM = dt_dL_dqdot - dL_dq - Q_ext;

% --- 3. Linearization ---
vars = [sym('v_dot'); sym('omega_dot')];
sol = solve(EOM == 0, vars);
f = [v; sol.v_dot; omega; sol.omega_dot];

states = [p; v; theta; omega];
inputs = [F1; F2]; % Note the two-column B matrix now

A_sym = jacobian(f, states);
B_sym = jacobian(f, inputs);

% Substitute parameters
params = {M, m, l, I, g, p, v, theta, omega};
values = {0.5, 0.2, 0.3, 0.006, 9.81, 0, 0, 0, 0}; 

A = double(subs(A_sym, params, values));
B = double(subs(B_sym, params, values));

% --- 4. Multi-Input LQR Design ---
Q = diag([100, 10, 100, 10]); 
% R must now be a 2x2 matrix because we have 2 inputs
R = diag([0.1, 0.1]); 

K = lqr(A, B, Q, R);

disp('State-Space B Matrix (Two Inputs):');
disp(B);
disp('LQR Gain Matrix K (2x4):');
disp(K);