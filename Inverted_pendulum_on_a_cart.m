clear all, clc

%% Initial Conditions

% Linearization point
q0 = [0; 0];
q_dot0 = [0; 0];

% Model conditions
tspan = 0:.001:10;
x0 = [-1; 0.5; 0; 0] + double([q0; q_dot0]);
wr = [1; 0; 0; 0] + double([q0; q_dot0]); % desired position

% Motors restrictions
F_max = 1000; % Max Newton or Nm your motor can provide

%% Define constants

% Floor
floor = struct;
floor.length = 5;
floor.width = 0.4;
floor.height = 0.01;

% Cart
cart = struct;
cart.length = 0.05;
cart.width = 0.05;
cart.height = 0.05;
cart.mass = 1;
cart.cof = 0; % coefficient of friction

% Rod
rod = struct;
rod.length = 0.3;
rod.width = 0.02;
rod.thickness = 0.005;
rod.mass = 1;
rod.cof = 0; % coefficient of friction

% Center of Mass
COM = struct;
COM.mass = rod.mass;

% Other parameters
motion_tc = 0.02;
g = 9.81;

total_mass = cart.mass + rod.mass

%% Define variables

% Cart
cart.x = sym('cart_x', 'real');
cart.x_dot = sym('cart_x_dot', 'real');
cart.x_ddot = sym('cart_x_ddot', 'real');

cart.F = sym('cart_F', 'real');

% Center of Mass
COM.theta = sym('COM_theta', 'real');
COM.theta_dot = sym('COM_theta_dot', 'real');
COM.theta_ddot = sym('COM_theta_ddot', 'real');

% States
q = [cart.x; COM.theta];
q_dot = [cart.x_dot; COM.theta_dot];
q_ddot = [cart.x_ddot; COM.theta_ddot];

%% Derive other parameters

% Cart
cart.y = cart.height;
cart.y_dot = jacobian(cart.y, q) * q_dot;

% Center of Mass
COM.l = rod.length / 2;
COM.l_dot = jacobian(COM.l, q) * q_dot;

COM.x = cart.x - COM.l * sin(COM.theta);
COM.x_dot = jacobian(COM.x, q) * q_dot;
COM.y = cart.y + COM.l * cos(COM.theta);
COM.y_dot = jacobian(COM.y, q) * q_dot;

% Inputs
u = [cart.F; 0];
u_max = [F_max; 0];

%% Define Lagrange Equations

% Compute kinetic energy
KE = 1/2 * cart.mass * (cart.x_dot ^ 2 + cart.y_dot ^ 2) + ...
     1/2 * COM.mass * (COM.x_dot ^ 2 + COM.y_dot ^ 2);

% Compute potential energy
PE = (COM.mass * COM.y + cart.mass * cart.y) * g;

%% Solve Lagrange Equations

L = KE - PE;

R = 1/2 * cart.cof * cart.x_dot ^ 2;
R = R + 1/2 * rod.cof * COM.theta_dot ^ 2;

% Compute the equations of motion using Lagrange's equations
EOM = jacobian(jacobian(L, q_dot), [q; q_dot]) * [q_dot; q_ddot] - jacobian(L, q)' + jacobian(R, q_dot)'

% Mass matrix D(q): coefficients of accelerations in EOM (linear in q_ddot)
D = jacobian(EOM, q_ddot);  % n x n

% Remaining terms: move D*q_ddot to left, Cg contains velocity and gravity and -u
Cg = EOM - D*q_ddot;        % n x 1

% Gravity vector G(q): set velocities to zero to isolate q-only terms
zeroDQ = sym(zeros(size(q_dot)));
Gvec = subs(Cg, q_dot, zeroDQ);  % n x 1

% Velocity-dependent terms C(q,q_dot) (Coriolis/centrifugal + other velocity terms)
Cvec = Cg - Gvec;  % n x 1

% Solve for nonlinear accelerations (q_ddot = D^{-1}*( -C + u ))
acc_nl = simplify(D \ (-Cg + u))   % n x 1 symbolic q_ddot expressions

%% Linearization
% Linearization about equilibrium (q0, q_dot0). Use symbolic q0,q_dot0 or numeric later.

% Evaluate D at equilibrium
D0 = subs(D, q, q0);  % D evaluated at q0 (no q_dot dependence)

% Linearize C: keep first-order in q_dot -> C_lin = (∂C/∂q_dot)|0 * q_dot
C_q_dot_jac = jacobian(Cvec, q_dot);   % n x n
C_lin = C_q_dot_jac * q_dot;           % n x 1 (first-order in q_dot)

% Linearize G: G_lin = (∂G/∂q)|0 * (q - q0)
G_q_jac = jacobian(Gvec, q);    % n x n
% use small displacement delta_q = q - q0; here we keep symbolic q and assume q0 zeros
delta_q = q - q0;
G_lin = G_q_jac * delta_q;      % n x 1 (first-order in q)

% Linearized implicit dynamics: D0*q_ddot + C_lin + G_lin = Q_lin
% For input linearization, linearize u if needed (here assume u linear in input F: u = u)
% If u depends on q or q_dot include jacobian terms similarly. For simplicity assume u = u:
% Solve for linear q_ddot: q_ddot = D0 \ ( -C_lin - G_lin + u )
q_ddot_lin = simplify(D0 \ ( -C_lin - G_lin + u ));  % n x 1 (affine in q, q_dot, F)

% Compute A,B matrices symbolically
A_lin_sym = simplify(jacobian([q_dot; q_ddot_lin], [q; q_dot]))   % 2n x 2n
B_lin_sym = simplify(jacobian([q_dot; q_ddot_lin], symvar(u))) % 2n x m, more robust below

% Evaluate A,B at equilibrium (substitute q->q0, q_dot->q_dot0)
A_lin = simplify(subs(A_lin_sym, [q; q_dot], [q0; q_dot0]))
B_lin = simplify(subs(B_lin_sym, [q; q_dot], [q0; q_dot0]))

A_lin = double(A_lin);
B_lin = double(B_lin);

%% Design LQR controller
Q = diag(10*ones(size([q; q_dot])));
R = diag(0.1*ones(size(symvar(u))));

% Q = diag([5 10 10 0.1]);
% R = diag([10]);

K = lqr(A_lin, B_lin, Q, R); % N = 0

disp('LQR Gain Matrix K:');
disp(K);

%% Simulate closed-loop system
u_law = @(x) max(-u_max, min(u_max, -K*(x - wr))); % control law

D_handle  = matlabFunction(D,  'vars', {q});
Cg_handle = matlabFunction(Cg, 'vars', {[q; q_dot]});

[t,x] = ode23tb(@(t, x) my_non_linear_model(t, x, u_law(x), D_handle, Cg_handle), tspan, x0);

figure
u_history = max(-u_max, min(u_max, -K*(x' - wr)));
plot(t, [x, u_history']);
legend('x', '\theta', 'v', '\omega', 'F', '\tau');

function [x_dot, u] = my_non_linear_model(t, x, u, D_func, Cg_func)
    q_i  = x(1:2);
    q_dot_i = x(3:4);

    D_val  = D_func(q_i); 
    Cg_val = Cg_func([q_i; q_dot_i]);

    q_ddot = D_val \ (u - Cg_val);

    x_dot = [q_dot_i;  q_ddot];
end