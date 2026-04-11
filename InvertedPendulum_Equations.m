clear all, close all, clc

%% Define constants
wheel = struct;
wheel.radius = 0.2 / 2;
wheel.thickness = 0.01;
wheel.mass = 10;
wheel.I = 1/2 * wheel.mass * wheel.radius ^ 2;

lowerLeg = struct;
lowerLeg.length = 0.3;
lowerLeg.width = 0.02;
lowerLeg.thickness = 0.005;

upperLeg = struct;
upperLeg.length = 0.3;
upperLeg.width = 0.02;
upperLeg.thickness = 0.005;

joint_offset = lowerLeg.length * 0.05;

COM = struct;
COM.mass = 10;

g = 9.81;

%% Lagrange Equations

% variables
syms t tau;
assume(t, 'real');

wheel.theta_ddot = sym('wheel_theta_ddot');
COM.l_ddot = sym('COM_l_ddot');

wheel.theta = symfun('wheel_theta(t)', t);
COM.l = symfun('COM_l(t)', t);
COM.theta = symfun('COM_theta(t)', t);

wheel.theta_dot = diff(wheel.theta, t);
wheel.x = wheel.theta * wheel.radius;
wheel.y = wheel.radius;
wheel.x_dot = wheel.theta_dot * wheel.radius;
wheel.y_dot = diff(wheel.y, t);

COM.l_dot = diff(COM.l, t);
COM.theta_dot = diff(COM.theta, t);
COM.x = wheel.x + COM.l * sin(COM.theta);
COM.y = wheel.radius + COM.l * cos(COM.theta);
COM.x_dot = diff(COM.x, t);
COM.y_dot = diff(COM.y, t);

% Compute kinetic energy
KE = 2 * (1/2 * wheel.mass * (wheel.x_dot ^ 2 + wheel.y_dot ^ 2)) + ...
     2 * (1/2 * wheel.I * wheel.theta_dot ^ 2) + ...
         (1/2 * COM.mass * (COM.x_dot ^ 2 + COM.y_dot ^ 2));

% Compute potential energy
PE = COM.mass * g * COM.y;

L = KE - PE;

q = [wheel.theta(t); COM.l(t)];
q_dot = [wheel.theta_dot(t); COM.l_dot(t)];
q_ddot = [wheel.theta_ddot; COM.l_ddot];
u = [tau; 0];

% Compute the equations of motion using Lagrange's equations
EOM = jacobian(jacobian(L, q_dot), [q; q_dot]) * [q_dot; q_ddot] - jacobian(L, q)' - u;

% EOM1 = [diff(diff(L, q_dot(1)), t) - diff(L, q(1)) - u(1);
%        diff(diff(L, q_dot(2)), t) - diff(L, q(2)) - u(2)];

q_ddot = struct2cell(solve(EOM == 0, q_ddot));
q_ddot = [q_ddot{:}].';

%% Linearization
states = [q; q_dot];
f = [q_dot; q_ddot];

% states = [q(1); q_dot(1); q(2); q_dot(2)];
% f = [q_dot(1); q_ddot(1); q_dot(2); q_ddot(2)];

A_sym = jacobian(f, states)
B_sym = jacobian(f, nonzeros(u));

% Initial conditions: wheel.theta(t), COM.l(t), then their first derivatives
x0 = [0; wheel.radius+lowerLeg.length-joint_offset; pi; 0];

A = subs(A_sym, states, x0);
val = 7.328745;
test_element = A(3,4);
simplify(subs(A, [COM.theta(t), conj(COM.theta(t))], 100*[pi/18, pi/18]))
vpa(A);


% A = subs(A_sym, [COM.theta_dot(t); states], [3; x0])
% B = subs(B_sym, [COM.theta_dot(t); states], [3; x0]);
% 
% A = double(A)
% B = double(B)
% A = [0 1 0 0;
%      0 -d/M b*m*g/M 0;
%      0 0 0 1;
%      0 -b*d/(M*L) -b*(m+M)*g/(M*L) 0];
% B = [0; 1/M; 0; b*1/(M*L)];

%% Design LQR controller
% Q = diag([100, 10, 100, 10]);
% R = diag([0.1]); 
% 
% K = lqr(A, B, Q, R); % N = 0
% 
% disp('LQR Gain Matrix K (2x4):');
% disp(K);

%% Simulate closed-loop system
% tspan = 0:.001:10;
% x0 = [-1; 0; pi+.1; 0]; % initial condition
% wr = [1; 0; pi; 0]; % reference position
% u=@(x)-K*(x - wr); % control law
% [t,x] = ode45(@(t,x)pendcart(x,m,M,L,g,d,u(x)),tspan,x0);
% 
% plot(t, x);
% legend('x', 'v', '\theta', '\omega');