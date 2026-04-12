clear all, close all, clc

%% Define constants

robot.width = 0.2;

% Floor
floor = struct;
floor.length = robot.width * 2;
floor.width = 5;
floor.thickness = 0.1;

% Wheel
wheel = struct;
wheel.radius = 0.2 / 2;
wheel.thickness = 0.01;
wheel.mass = 1;
wheel.I = 1/2 * wheel.mass * wheel.radius ^ 2;

% Upper Leg
lowerLeg = struct;
lowerLeg.length = 0.3;
lowerLeg.width = 0.02;
lowerLeg.thickness = wheel.thickness * 1.5; % 0.005
lowerLeg.mass = 1;

% Upper Leg
upperLeg = struct;
upperLeg.length = 0.3;
upperLeg.width = 0.02;
upperLeg.thickness = 0.005;
upperLeg.mass = 0; % 1

% Center of Mass
COM = struct;
COM.mass = lowerLeg.mass;

% Other parameters
motion_tc = 0.02;
joint_offset = lowerLeg.length * 0.05;
g = 9.81;
payload = 0; % 15

total_mass = 2 * (wheel.mass + lowerLeg.mass + upperLeg.mass) + payload


%% Lagrange Equations

% variables
syms tau real

wheel.theta = sym('wheel_theta', 'real');
wheel.theta_dot = sym('wheel_theta_dot', 'real');
wheel.theta_ddot = sym('wheel_theta_ddot', 'real');

wheel.x = wheel.theta * wheel.radius;
wheel.y = wheel.radius;
wheel.x_dot = wheel.theta_dot * wheel.radius;
wheel.y_dot = 0;

COM.l = wheel.radius+lowerLeg.length-joint_offset;
COM.l_dot = 0;
% COM.l = sym('COM_l', 'real');
% COM.l_dot = sym('COM_l_dot', 'real');
% COM.l_ddot = sym('COM_l_ddot', 'real');
COM.theta = sym('COM_theta', 'real');
COM.theta_dot = sym('COM_theta_dot', 'real');
COM.theta_ddot = sym('COM_theta_ddot', 'real');

COM.x = wheel.x + COM.l * sin(COM.theta);
COM.y = wheel.radius + COM.l * cos(COM.theta);
COM.x_dot = wheel.x_dot + COM.l * cos(COM.theta) * COM.theta_dot + sin(COM.theta) * COM.l_dot;
COM.y_dot = - COM.l * sin(COM.theta) * COM.theta_dot + cos(COM.theta) * COM.l_dot;

% Compute kinetic energy
KE = 2 * (1/2 * wheel.mass * (wheel.x_dot ^ 2 + wheel.y_dot ^ 2)) + ...
     2 * (1/2 * wheel.I * wheel.theta_dot ^ 2) + ...
         (1/2 * COM.mass * (COM.x_dot ^ 2 + COM.y_dot ^ 2));

% Compute potential energy
PE = COM.mass * g * COM.y;

L = KE - PE;

q = [wheel.theta; COM.theta];
q_dot = [wheel.theta_dot; COM.theta_dot];
q_ddot = [wheel.theta_ddot; COM.theta_ddot];
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

A_sym = simplify(jacobian(f, states));
B_sym = simplify(jacobian(f, nonzeros(u)));

% Initial conditions: wheel.theta, COM.l, COM.theta, then their first derivatives
% x0 = [0; wheel.radius+lowerLeg.length-joint_offset; pi/18; pi; 0; 0];
% Initial conditions: wheel.theta, COM.theta, then their first derivatives
% x0 = [0; 0; 0; 0];
x0 = [0; 0; pi/2; pi/18];

A = simplify(subs(A_sym, states, x0));
B = subs(B_sym, states, x0);

% A = subs(A, [COM.theta; COM.theta_dot], [pi/18; 0]);
% B = subs(B, [COM.theta; COM.theta_dot], [pi/18; 0]);

B = double(B);
A = double(A);

%% Design LQR controller
Q = diag([100, 10, 100, 10]);
R = diag([0.1]); 

K = lqr(A, B, Q, R); % N = 0

disp('LQR Gain Matrix K:');
disp(K);

%% Simulate closed-loop system
% tspan = 0:.001:10;
% x0 = [-1; 0; pi+.1; 0]; % initial condition
wr = [wheel.radius+floor.length/3; 0; 0; 0]; % reference position
% u=@(x)-K*(x - wr); % control law
% [t,x] = ode45(@(t,x)pendcart(x,m,M,L,g,d,u(x)),tspan,x0);
% 
% plot(t, x);
% legend('x', 'v', '\theta', '\omega');













% (1962*cos(pi/18)^2 - 1962*sin(pi/18)^2)/(20*(3*cos(pi/18)^2 + 4*sin(pi/18)^2)) - (cos(pi/18)*sin(pi/18)*(200*tau*cos(pi/18)^2 + 200*tau*sin(pi/18)^2 - 1962*cos(pi/18)*sin(pi/18)))/(10*(3*cos(pi/18)^2 + 4*sin(pi/18)^2)^2)
% (7848*cos(pi/18) + 200*tau*sin(pi/18))/(77*(3*cos(pi/18)^2 + 4*sin(pi/18)^2)) - (2*cos(pi/18)*sin(pi/18)*(7848*sin(pi/18) - 200*tau*cos(pi/18)))/(77*(3*cos(pi/18)^2 + 4*sin(pi/18)^2)^2)




