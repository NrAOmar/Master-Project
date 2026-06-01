% clear all, clc

%% Initial Conditions

tic
% Linearization point
q0 = [0; 0];
q_dot0 = [0; 0];

% Model conditions
tspan = 0:.001:30;
wr = [2/0.1; 0; 0; 0] + double([q0; q_dot0]); % desired position
x0 = [-wr(1); 0; 0; 0] + double([q0; q_dot0]);

% Motors restrictions
tau_max = 20; % Max Newton or Nm your motor can provide

required_height = 0.8;

%% Define constants

% Floor
floor = struct;
floor.length = 10;
floor.width = 0.4;
floor.height = 0.01;

% Wheel
wheel = struct;
wheel.radius = 0.2 / 2;
wheel.thickness = 0.1 * wheel.radius;
wheel.mass = 0.5;
wheel.I = 1/2 * wheel.mass * (wheel.radius ^ 2 + (wheel.radius - wheel.thickness)^ 2);

% Rod
rod = struct;
rod.length = (required_height - wheel.radius) / 2;
rod.width = 0.02;
rod.thickness = 0.005;
rod.mass = 3;

lowerLeg = rod;
upperLeg = rod;

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

% Lower Leg
lowerLeg.theta = sym('lowerLeg_theta', 'real');
lowerLeg.theta_dot = sym('lowerLeg_theta_dot', 'real');
lowerLeg.theta_ddot = sym('lowerLeg_theta_ddot', 'real');

% States
q = [wheel.theta; lowerLeg.theta];
q_dot = [wheel.theta_dot; lowerLeg.theta_dot];
q_ddot = [wheel.theta_ddot; lowerLeg.theta_ddot];

%% Derive other parameters

% Wheel
wheel.x = wheel.radius * wheel.theta;
wheel.x_dot = jacobian(wheel.x, q) * q_dot;
wheel.y = wheel.radius;
wheel.y_dot = jacobian(wheel.y, q) * q_dot;

% Lower Leg
lowerLeg.speed_coef = 0; % coefficient of friction
lowerLeg.length = (required_height - wheel.radius) * 0.3;
lowerLeg.x = wheel.x + lowerLeg.length / 2 * cos(pi/2 + lowerLeg.theta);
lowerLeg.x_dot = jacobian(lowerLeg.x, q) * q_dot;
lowerLeg.y = wheel.y + lowerLeg.length / 2 * sin(pi/2 + lowerLeg.theta);
lowerLeg.y_dot = jacobian(lowerLeg.y, q) * q_dot;

% Upper Leg
upperLeg.speed_coef = 0; % coefficient of friction
upperLeg.length = (required_height - wheel.radius) * 0.7;
upperLeg.theta = sym('upperLeg_theta', 'real');
upperLeg.x = lowerLeg.x + lowerLeg.length / 2 * cos(pi/2 + lowerLeg.theta) + upperLeg.length / 2 * cos(upperLeg.theta + lowerLeg.theta);
upperLeg.x_dot = jacobian(upperLeg.x, q) * q_dot;
upperLeg.y = lowerLeg.y + lowerLeg.length / 2 * sin(pi/2 + lowerLeg.theta) + upperLeg.length / 2 * sin(upperLeg.theta + lowerLeg.theta);
upperLeg.y_dot = jacobian(upperLeg.y, q) * q_dot;
upperLeg.theta0 = 0;

% Payload
payload.x = upperLeg.x + upperLeg.length / 2 * cos(upperLeg.theta + lowerLeg.theta);
payload.x_dot = jacobian(payload.x, q) * q_dot;
payload.y = upperLeg.y + upperLeg.length / 2 * sin(upperLeg.theta + lowerLeg.theta);
payload.y_dot = jacobian(payload.y, q) * q_dot;

% Center of mass
COM.x = (2 * lowerLeg.x * lowerLeg.mass + 2 * upperLeg.x * upperLeg.mass + payload.x * payload.mass) / COM.mass;
COM.x_dot = jacobian(COM.x, q) * q_dot;
COM.y = (2 * lowerLeg.y * lowerLeg.mass + 2 * upperLeg.y * upperLeg.mass + payload.y * payload.mass) / COM.mass;
COM.y_dot = jacobian(COM.y, q) * q_dot;
COM.l = simplify(sqrt((COM.x - wheel.x) ^ 2 + (COM.y - wheel.y) ^ 2));
COM.theta = atan2((COM.y - wheel.y), (COM.x - wheel.x)) - pi/2;

matlabFunctionBlock( ...
    "Balancing_Robot_model/Balancing Robot/Display COM/calculate_COM_x", ...
    COM.x, ...
    'Vars',[wheel.theta, upperLeg.theta, lowerLeg.theta]);

matlabFunctionBlock( ...
    "Balancing_Robot_model/Balancing Robot/Display COM/calculate_COM_y", ...
    COM.y, ...
    'Vars',[wheel.theta, upperLeg.theta, lowerLeg.theta]);

matlabFunctionBlock( ...
    "Balancing_Robot_model/Sensor/calculate_COM_theta", ...
    subs(COM.theta, 'lowerLeg_theta', q0(2)), ...
    'Vars',[upperLeg.theta]);

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

R =     1/2 * lowerLeg.speed_coef * wheel.theta_dot ^ 2;
R = R + 1/2 * upperLeg.speed_coef * lowerLeg.theta_dot ^ 2;

% Compute the equations of motion using Lagrange's equations
EOM = jacobian(jacobian(L, q_dot), [q; q_dot]) * [q_dot; q_ddot] - jacobian(L, q)' + jacobian(R, q_dot)';

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
acc_nl = simplify(D \ (-Cg + u));   % n x 1 symbolic q_ddot expressions

toc
tic

K = zeros(1, 4, round((pi/180 / pi/2)) + 1);
for i = 0:90

upperLeg.theta0 = i * pi / 180;
q0(2) = 0;

%% Linearization
% Linearization about equilibrium (q0, q_dot0). Use symbolic q0,q_dot0 or numeric later.
% Balance at center of mass angle instead of payload angle

q0(2) = q0(2) - double(subs(COM.theta, {'lowerLeg_theta', 'upperLeg_theta'}, {q0(2), upperLeg.theta0}));
x0 = x0 + double([q0; q_dot0]);

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
A_lin_sym = simplify(jacobian([q_dot; q_ddot_lin], [q; q_dot]));   % 2n x 2n
B_lin_sym = simplify(jacobian([q_dot; q_ddot_lin], symvar(u))); % 2n x m, more robust below

% Evaluate A,B at equilibrium (substitute q->q0, q_dot->q_dot0)
A_lin = simplify(subs(A_lin_sym, [q; q_dot], [q0; q_dot0]));
B_lin = simplify(subs(B_lin_sym, [q; q_dot], [q0; q_dot0]));

A_lin = double(subs(A_lin, 'upperLeg_theta', 0));
B_lin = double(subs(B_lin, 'upperLeg_theta', 0));

toc
tic

%% Design LQR controller
Q = diag(((1 ./ [wr(1)*0.05 0.01 wr(1)*0.2 0.2]) .^ 2));
R = diag(((1 ./ nonzeros(u_max)) .^ 2));

K(:,:,i+1) = lqr(A_lin, B_lin, Q, R); % N = 0

disp('For Upper Leg theta:');
disp(upperLeg.theta0 * 180 / pi)
disp('LQR Gain Matrix K:');
disp(K);

end
toc

%% Simulate closed-loop system

% u_law = @(x) max(-u_max, min(u_max, -K*(x - wr))); % control law
% 
% D_handle  = matlabFunction(D,  'vars', {q});
% Cg_handle = matlabFunction(Cg, 'vars', {[q; q_dot]});
% 
% [t,x] = ode23tb(@(t, x) my_non_linear_model(t, x, u_law(x), D_handle, Cg_handle), tspan, x0);
% 
% figure
% u_history = max(-u_max, min(u_max, -K*(x' - wr)));
% plot(t, [x]);
% legend('\theta_{wheel}', '\theta_{rod}', '\omega_{wheel}', '\omega_{rod}', '\tau_{wheel}', '\tau_{rod}');
% 
% function [x_dot, u] = my_non_linear_model(t, x, u, D_func, Cg_func)
%     q_i  = x(1:numel(x)/2);
%     q_dot_i = x(numel(x)/2+1:end);
% 
%     D_val  = D_func(q_i); 
%     Cg_val = Cg_func([q_i; q_dot_i]);
% 
%     q_ddot = D_val \ (u - Cg_val);
% 
%     x_dot = [q_dot_i;  q_ddot];
% end

%% Height Control

syms theta(t) T(t) theta_s T_s s

upperLeg.theta_dot = diff(theta,t,1);
upperLeg.theta_ddot = diff(theta,t,2);

upperLeg.theta_offset = pi/2;
theta_lower = -subs(COM.theta, {'lowerLeg_theta', 'upperLeg_theta'}, {0, theta});

% 2. Linearize @ theta = 0
% ode = subs(ode, sin(theta), theta);
% ode = subs(ode, cos(theta), 1);

theta0_val = 0;
f = cos(theta + theta_lower);

% compute zeroth and first derivatives
f0 = simplify(subs(f, theta, theta0_val));
df = simplify(diff(f, theta));
df0 = simplify(subs(df, theta, theta0_val));

% linearized expression: f ≈ f0 + df0*(theta - theta0)
cos_lin = simplify(f0 + df0*theta);

% 1. Physical ODE
% upperLeg.Torque = ...
%       payload.mass * upperLeg.length * cos(theta + theta_lower) * g / 2 + ... % for 2 joints
%       upperLeg.mass * (upperLeg.length / 2) * cos(theta + theta_lower) * g + ...
%       upperLeg.mass * (upperLeg.length / 2) ^ 2 * upperLeg.theta_ddot + ...
%       upperLeg.speed_coef * upperLeg.theta_dot;

upperLeg.Torque = ...
      payload.mass * upperLeg.length * cos_lin * g / 2 + ... % for 2 joints
      upperLeg.mass * (upperLeg.length / 2) * cos_lin * g + ...
      upperLeg.mass * (upperLeg.length / 2) ^ 2 * upperLeg.theta_ddot + ...
      upperLeg.speed_coef * upperLeg.theta_dot;

upperLeg.K = abs(double(subs(upperLeg.Torque, {theta, upperLeg.theta_dot, upperLeg.theta_ddot}, {0, 0, 0}) / upperLeg.theta_offset));
upperLeg.Torque = upperLeg.Torque - upperLeg.K * (upperLeg.theta_offset - theta);

ode = upperLeg.Torque == T(t);

% 3. Laplace and clear initial conditions
ode_s = simplify(laplace(ode));
L_ode = subs(ode_s, {laplace(theta(t)), laplace(T(t)), theta(0), subs(diff(theta,t),t,0)}, {theta_s, T_s, 0, 0});

% 4. Solve for Xs/Fs
G = simplify(solve(L_ode, theta_s) / T_s);

% get coefficient vectors in descending powers of s
[num_sym, den_sym] = numden(G);
num_coeffs = fliplr(coeffs(num_sym, s, 'All'));   % descending order
den_coeffs = fliplr(coeffs(den_sym, s, 'All'));

% create MATLAB functions that return numeric coefficient vectors for a given T_s
num_fun = matlabFunction(num_coeffs, 'Vars', T_s);
den_fun = matlabFunction(den_coeffs, 'Vars', T_s);

% when you pick a value for T_s (e.g. Tval)
Tval = 0.01;                      % choose a candidate T_s
num = double(num_fun(Tval));
den = double(den_fun(Tval));
plant = tf(num, den);             % numeric transfer function you can feed to pidtune

% Create a tuning option for a fast response
opts = pidtuneOptions( ...
    'PhaseMargin', 65, ...         % Lower margin = faster, but more overshoot
    'DesignFocus', 'disturbance-rejection'); % Focus on staying still when pushed

% [C_pid, info] = pidtune(plant, 'pid', opts);
[C_pid, info] = pidtune(plant, 'pid');
% C_pid = pidTuner(plant, 'pid');

%% See how it performs

% sys_cl = feedback(C_pid * plant, 1);
% figure;
% step(sys_cl)
% title('PID Control Response')
% grid on