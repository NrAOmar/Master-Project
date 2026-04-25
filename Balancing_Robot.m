clear all, clc

%% Initial Conditions

% Linearization point
q0 = [0; 0];
q_dot0 = [0; 0];

% Model conditions
tspan = 0:.001:10;
x0 = [0; 0; 0; 0] + double([q0; q_dot0]);
wr = [6*pi; 0; 0; 0] + double([q0; q_dot0]); % desired position

% Motors restrictions
tau_max = 1000; % Max Newton or Nm your motor can provide

required_height = 0.8;

%% Define constants

% Floor
floor = struct;
floor.length = 5;
floor.width = 0.4;
floor.height = 0.01;

% Wheel
wheel = struct;
wheel.radius = 0.2 / 2;
wheel.thickness = 0.1 * wheel.radius;
wheel.mass = 0.5;
wheel.I = 1/2 * wheel.mass * (wheel.radius ^ 2 + (wheel.radius - wheel.thickness)^ 2);
wheel.cof = 0; % coefficient of friction

% Rod
rod = struct;
rod.length = (required_height - wheel.radius) / 2;
rod.width = 0.02;
rod.thickness = 0.005;
rod.mass = 3;
rod.cof = 0; % coefficient of friction

% Payload
payload = struct;
payload.mass = 10;

% Center of Mass
COM = struct;
COM.mass = 4 * rod.mass + payload.mass;

% Other parameters
motion_tc = 0.02;
g = 9.80665;

total_mass = 2 * wheel.mass + COM.mass;

%% Define variables

% Wheel
wheel.theta = sym('wheel_theta', 'real');
wheel.theta_dot = sym('wheel_theta_dot', 'real');
wheel.theta_ddot = sym('wheel_theta_ddot', 'real');

wheel.tau = sym('wheel_tau', 'real');

% Center of Mass
COM.theta = sym('COM_theta', 'real');
COM.theta_dot = sym('COM_theta_dot', 'real');
COM.theta_ddot = sym('COM_theta_ddot', 'real');

% States
q = [wheel.theta; COM.theta];
q_dot = [wheel.theta_dot; COM.theta_dot];
q_ddot = [wheel.theta_ddot; COM.theta_ddot];

%% Derive other parameters

% Wheel
wheel.x = wheel.radius * wheel.theta;
wheel.x_dot = jacobian(wheel.x, q) * q_dot;
wheel.y = wheel.radius;
wheel.y_dot = jacobian(wheel.y, q) * q_dot;

% Lower Leg
lowerLeg = rod;
lowerLeg.theta = pi/6;
lowerLeg.theta_dot = jacobian(lowerLeg.theta, q) * q_dot;
lowerLeg.x = wheel.x + lowerLeg.length / 2 * cos(pi/2 + lowerLeg.theta + COM.theta);
lowerLeg.x_dot = jacobian(lowerLeg.x, q) * q_dot;
lowerLeg.y = wheel.y + lowerLeg.length / 2 * sin(pi/2 + lowerLeg.theta + COM.theta);
lowerLeg.y_dot = jacobian(lowerLeg.y, q) * q_dot;

% Upper Leg
upperLeg = rod;
upperLeg.theta = -lowerLeg.theta;
upperLeg.theta_dot = jacobian(upperLeg.theta, q) * q_dot;
upperLeg.x = lowerLeg.x + lowerLeg.length / 2 * cos(pi/2 + lowerLeg.theta + COM.theta) + upperLeg.length / 2 * cos(pi/2 + upperLeg.theta + COM.theta);
upperLeg.x_dot = jacobian(upperLeg.x, q) * q_dot;
upperLeg.y = lowerLeg.y + lowerLeg.length / 2 * sin(pi/2 + lowerLeg.theta + COM.theta) + upperLeg.length / 2 * sin(pi/2 + upperLeg.theta + COM.theta);
upperLeg.y_dot = jacobian(upperLeg.y, q) * q_dot;

% Payload
payload.l = lowerLeg.length * sin(pi/2 + lowerLeg.theta) + upperLeg.length * sin(pi/2 + upperLeg.theta);
payload.l_dot = jacobian(payload.l, q) * q_dot;
payload.x = upperLeg.x + upperLeg.length / 2 * cos(pi/2 + upperLeg.theta + COM.theta);
payload.x_dot = jacobian(payload.x, q) * q_dot;
payload.y = upperLeg.y + upperLeg.length / 2 * sin(pi/2 + upperLeg.theta + COM.theta);
payload.y_dot = jacobian(payload.y, q) * q_dot;

% Center of mass
COM.x = (2 * lowerLeg.x * lowerLeg.mass + 2 * upperLeg.x * upperLeg.mass + payload.x * payload.mass) / COM.mass;
COM.x_dot = jacobian(COM.x, q) * q_dot;
COM.y = (2 * lowerLeg.y * lowerLeg.mass + 2 * upperLeg.y * upperLeg.mass + payload.y * payload.mass) / COM.mass;
COM.y_dot = jacobian(COM.y, q) * q_dot;

% Balance at center of mass angle instead of payload angle
COM.theta0 = double(subs(atan2(-(COM.x-wheel.x), (COM.y-wheel.y)), 'COM_theta', q0(2)));
q0(2) = q0(2) - COM.theta0;
x0 = x0 + double([q0; q_dot0]);
wr = wr + double([q0; q_dot0]);

% Inputs
u = [wheel.tau; 0];
u_max = [tau_max; 0];

%% Define Lagrange Equations

% Compute kinetic energy
KE = 2 * (1/2 * wheel.mass * (wheel.x_dot ^ 2 + wheel.y_dot ^ 2)) + ...
     2 * (1/2 * wheel.I * wheel.theta_dot ^ 2) + ...
          1/2 * COM.mass * (COM.x_dot ^ 2 + COM.y_dot ^ 2);

% Compute potential energy
PE = (COM.mass * COM.y + 2 * wheel.mass * wheel.y) * g;

%% Solve Lagrange Equations

L = KE - PE;

R =     1/2 * wheel.cof * wheel.theta_dot ^ 2;
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

% Q = diag([1 100 100 100]);
% R = diag([1]);

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
plot(t, [x]);
legend('\theta_{wheel}', '\theta_{rod}', '\omega_{wheel}', '\omega_{rod}', '\tau_{wheel}', '\tau_{rod}');

function [x_dot, u] = my_non_linear_model(t, x, u, D_func, Cg_func)
    q_i  = x(1:numel(x)/2);
    q_dot_i = x(numel(x)/2+1:end);

    D_val  = D_func(q_i); 
    Cg_val = Cg_func([q_i; q_dot_i]);

    q_ddot = D_val \ (u - Cg_val);

    x_dot = [q_dot_i;  q_ddot];
end