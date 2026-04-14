function sf_aerodyn(block)
% S-function sf_aerodyn.m
% This S-function represents the nonlinear aircraft dynamics

% Copyright 1986-2022 The MathWorks, Inc. 

    setup(block);
end

function setup(block)
    mdlInitializeSizes(block);
    block.RegBlockMethod('InitializeConditions',@mdlInitConditions);
    block.RegBlockMethod('Derivatives', @mdlDerivatives);
    block.RegBlockMethod('Outputs', @mdlOutputs);
end

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function mdlInitializeSizes(block)
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%

    % Register number of ports
    block.NumInputPorts = 1;
    block.NumOutputPorts = 1;
    
    % Setup port properties to be inherited or dynamic
    block.SetPreCompInpPortInfoToDynamic;
    block.SetPreCompOutPortInfoToDynamic;
    
    % Override port properties
    block.InputPort(1).DatatypeID = 0;
    block.InputPort(1).Dimensions = 4;
    block.InputPort(1).DirectFeedthrough = true;
    block.OutputPort(1).DatatypeID = 0;
    block.OutputPort(1).Dimensions = 8;

    % Register number of states
    block.NumContStates = 8;

    block.SimStateCompliance = 'DefaultSimState';
    
    % initialize the array of sample times
    block.SampleTimes = [0 0];
end

%
%=============================================================================
% mdlInitConditions
% Initialize continuous states.
%=============================================================================
%
function mdlInitConditions(block)
    % Initial conditions
    % x = [u,v,w,p,q,r,theta,phi]';
    % To linearized the model assume  phi = 0 and theta = 0; 
    block.ContStates.Data = [0 0 0 0 0 0 0 0];
end

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function mdlDerivatives(block)
    x = block.ContStates.Data;
    u = block.InputPort(1).Data;

    a1 = [-0.0404    0.0618    0.0501   -0.0000   -0.0005    0.0000 
       -0.1686   -1.1889    7.6870         0    0.0041         0
        0.1633   -2.6139   -3.8519    0.0000    0.0489   -0.0000 
       -0.0000   -0.0000   -0.0000   -0.3386   -0.0474   -6.5405  
       -0.0000    0.0000   -0.0000   -1.1288   -0.9149   -0.3679  
       -0.0000   -0.0000   -0.0000    0.9931   -0.1763   -1.2047 
             0         0    0.9056         0         0   -0.0000
             0         0   -0.0000         0    0.9467   -0.0046];
    a = [a1, zeros(8,2)];
     
    b1 =[ -0.0404    0.0618    0.0501
       -0.1686   -1.1889    7.6870
        0.1633   -2.6139   -3.8519
       -0.0000   -0.0000   -0.0000
       -0.0000    0.0000   -0.0000
       -0.0000   -0.0000   -0.0000
             0         0    0.9056
             0         0   -0.0000];
     
    b2 =[ 20.3929   -0.4694   -0.2392   -0.7126
        0.1269   -2.6932    0.0013    0.0033
      -64.6939  -75.6295    0.6007    3.2358
       -0.0000         0    0.1865    3.6625
       -0.0000         0   23.6053    5.6270
       -0.0001         0    3.9462  -41.4112
             0         0         0         0
             0         0         0         0];
    
    g=32.2;
    %% The state vector is x = [u,v,w,p,q,r,theta,phi]'%%  
    
    n1 = -g*sin(x(7));
    n2 = g*cos(x(7))*sin(x(8));
    n3 = g*cos(x(7))*cos(x(8));
    n7 = x(5)*cos(x(8))-x(6)*sin(x(8));
    n8 = (x(5)*sin(x(8)) + x(6)*cos(x(8)))*tan(x(7));
    
    % dx/dt = a*x + b1*w + b2*u + [n1;n2;n3;0;0;0;n7;n8];
    block.Derivatives.Data = a*x + b2*u + [n1;n2;n3;0;0;0;n7;n8];
end

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function mdlOutputs(block)
    block.OutputPort(1).Data = block.ContStates.Data;
    mdlGetTimeOfNextVarHit(block);
end

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function mdlGetTimeOfNextVarHit(block)
    sampleTime = 1;    %  Example, set the next hit to be one second later.
    block.NextTimeHit = block.CurrentTime + sampleTime;
end
