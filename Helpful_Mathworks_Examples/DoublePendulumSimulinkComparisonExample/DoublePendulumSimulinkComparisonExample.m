%% Double Pendulum in Simulink and Simscape Multibody
% 
% This example shows two models of a double pendulum, one using Simulink(R)
% input/output blocks and one using Simscape(TM) Multibody(TM).
% 
% The Simulink model uses signal connections, which define how data flows
% from one block to another.  The Simscape Multibody model is built using
% physical connections, which permit a bidirectional flow of energy between
% components.  Physical connections make it possible to add further stages
% to the pendulum simply by using copy and paste. Input/output connections
% require rederiving and reimplementing the equations.
% 
% The initial angle for each joint is defined by a MATLAB(R) variable.  The
% annotations on the Integrator blocks show the initial angles of the
% joints with respect to the world frame.
% 
% Copyright 2023 The MathWorks, Inc.

%% Model

mdl = 'DoublePendulumSimulinkComparison';
open_system(mdl)
set_param(mdl,'SimMechanicsOpenEditorOnUpdate','off');

%% Simulation Results from Simscape Logging

fh = PlotJointAngles;

%%
close(fh);
close_system(mdl, 0)