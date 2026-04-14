%% Bumper Car Parameters Initialization

% Copyright 2019-2023 The MathWorks, Inc.

%% gravity %%

g = 9.8; % m/s^2

%% Ramp structure %%

% Dimensions

ramp_structure.width        = 19;  % cm
ramp_structure.length       = 100; % cm
ramp_structure.ramp.angle   = 8;   % deg

% Graphic

ramp_structure.ramp.color   = [0.4196 0.5569 0.1373];
ramp_structure.rail.color   = [0.3961 0.3020 0.1647];
ramp_structure.pillar.color = [0.7569 0.6039 0.4196];
ramp_structure.base.color   = [0.2549 0.4078 0.1451];

%% Car %%

% Dimensions

car.body.length     = 14;  % cm
car.body.width      = 8;   % cm
car.body.height     = 3;   % cm
car.body.mass       = 0.4; % kg

car.wheel.radius    = 3;   % cm
car.wheel.thickness = 1.5; % cm
car.wheel.mass      = 0.1; % kg

car.start_dist      = 15;  % cm
car.start_height    = 3;   % cm
car.start_pitch     = 0;   % deg

% Graphic

car.body.color        = [0.3843 0.3216 0.2588];
car.wheel.color       = [0.6392 0.2353 0.0196];
car.wheel.spoke.color = [0.9608 0.8902 0.7333];

%% Contact force parameters %%

stiffness               = 10000; % N/m
damping                 = 30;    % N/(m/s)
transition_region_width = 1e-5;  % m
static_friction_coef    = 0.2;
kinetic_friction_coef   = 0.2;
critical_velocity       = 1e-2;  % m/s
