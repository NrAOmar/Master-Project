%% Computing Actuator Torques Using Inverse Dynamics
% This example illustrates the use of motion actuation to determine the
% actuator torques needed for the robot to achieve a given welding task.
% The system consists of a seven degree of freedom robot carrying a welding
% torch. The tip of the torch needs to trace the joints being welded. In
% this example the tip of the torch is made to trace (using motion
% actuation) a plus sign, a circle and a star sign on the workpiece. The
% torch is lifted off the workpiece when transitioning between the
% different shapes. The motion of the welding torch is specified and the
% actuator torques required at the various joints of the robot to achieve
% this motion is computed.
%
% Copyright 2023 The MathWorks, Inc.

open_system('InverseDynamics');
%%
close_system('InverseDynamics');