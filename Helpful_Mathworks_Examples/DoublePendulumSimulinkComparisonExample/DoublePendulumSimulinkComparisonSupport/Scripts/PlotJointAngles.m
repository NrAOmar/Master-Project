function fh = PlotJointAngles()
% Code to plot the joint angles from DoublePendulumSimulinkComparison

% Copyright 2023 The MathWorks, Inc.

sim('DoublePendulumSimulinkComparison');

% Get simulation results
Rz1q = simlogDoublePendulumSimulinkComparison.Revolute_Joint_1.Rz.q.series;
Rz2q = simlogDoublePendulumSimulinkComparison.Revolute_Joint_2.Rz.q.series;
theta1_sl = get(logsoutDoublePendulumSimulinkComparison,'theta1_sl');
theta2_in_world_sl = get(logsoutDoublePendulumSimulinkComparison,'theta2_in_world_sl');

% Plot simulation results
fh = figure('Name','DoublePendulumSimulinkComparison');
plot(Rz1q.time,Rz1q.values,'LineWidth',3)
hold on
plot(Rz2q.time,Rz2q.values,'LineWidth',3)
plot(theta1_sl.Values.time,theta1_sl.Values.Data*180/pi,'c--','LineWidth',1)
plot(theta1_sl.Values.time,(theta2_in_world_sl.Values.Data-theta1_sl.Values.Data)*180/pi,'c--','LineWidth',1)
hold off
grid on
title('Joint Angles');
ylabel('Angle (deg)');
xlabel('Time (s)');
legend({'Simscape Multibody Joint 1','Simscape Multibody Joint 2','Simulink'},'Location','Best');
