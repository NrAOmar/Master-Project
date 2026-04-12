clear all, clc

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

% Upper Leg
rod = struct;
rod.length = 0.3;
rod.width = 0.02;
rod.thickness = 0.005;
rod.mass = 1;

% Center of Mass
COM = struct;
COM.mass = rod.mass;

% Other parameters
motion_tc = 0.02;
g = 9.81;

total_mass = rod.mass

%% Define variables

% Cart
cart.x = sym('cart_x', 'real');
cart.x_dot = sym('cart_x_dot', 'real');
cart.x_ddot = sym('cart_x_ddot', 'real');
cart.y = cart.height;
cart.y_dot = 0;
cart.y_ddot = 0;

% Center of Mass
COM.theta = sym('COM_theta', 'real');
COM.theta_dot = sym('COM_theta_dot', 'real');
COM.theta_ddot = sym('COM_theta_ddot', 'real');
COM.l = rod.length / 2;
COM.l_dot = 0;
COM.l_ddot = 0;

COM.x = cart.x + COM.l * sin(COM.theta);
COM.x_dot = jacobian(COM.x, [cart.x; COM.theta]) * [cart.x_dot; COM.theta_dot];
COM.y = cart.y + COM.l * cos(COM.theta);
COM.y_dot = jacobian(COM.y, [COM.theta]) * [COM.theta_dot];

q = [cart.x; COM.theta];
q_dot = [cart.x_dot; COM.theta_dot];
q_ddot = [cart.x_ddot; COM.theta_ddot];

syms F real
u = [F; 0];

%% Define Lagrange Equations

% Compute kinetic energy
KE = 1/2 * cart.mass * (cart.x_dot ^ 2 + cart.y_dot ^ 2) + ...
     1/2 * COM.mass * (COM.x_dot ^ 2 + COM.y_dot ^ 2);

% Compute potential energy
PE = COM.mass * g * COM.y;

%% Solve Lagrange Equations

L = KE - PE;

b = 1;
R = 1/2 * b * cart.x_dot ^ 2;

% Compute the equations of motion using Lagrange's equations
EOM = jacobian(jacobian(L, q_dot), [q; q_dot]) * [q_dot; q_ddot] - jacobian(L, q)' + jacobian(R, q_dot)'

q_ddot_sol = struct2cell(solve(EOM == u, q_ddot));
q_ddot_sol = simplify([q_ddot_sol{:}].')

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
acc_nl = simplify(D \ (-Cvec + u))   % n x 1 symbolic q_ddot expressions

%% Linearization
% Linearization about equilibrium (q0, q_dot0). Use symbolic q0,q_dot0 or numeric later.
% Define equilibrium point here (example upright at zero)
q0 = sym(zeros(size(q)));    % change if needed
q_dot0 = sym(zeros(size(q_dot)));

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
% Q = diag([100, 10, 100, 10]);
% R = diag([0.1]); 
Q = eye(4); % 4x4 identify matrix
R = .0001;

K = lqr(A_lin, B_lin, Q, R); % N = 0

disp('LQR Gain Matrix K:');
disp(K);

%% Simulate closed-loop system
tspan = 0:.01:0.5;
% x0 = [-1; 0; pi+.1; 0]; % initial condition
x0 = [-1; -0.1; 0; 0];
wr = [floor.length/3; 0; 0; 0]; % reference position
u_controlled=@(x)-K*(x - wr); % control law

addpath(genpath('./From_Literature'))
[t,x] = ode23tb(@(t,x)pendcart(x,COM.mass,cart.mass,COM.l,g,0.5,u_controlled(x)),tspan,x0);

figure
plot(t, x);
legend('x', '\theta', 'v', '\omega');
ylim([-10 10]);