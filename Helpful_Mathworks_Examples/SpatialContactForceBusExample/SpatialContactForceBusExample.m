%% Using the Spatial Contact Force Block - Bumper Car
% This example shows a toy bumper car traveling down a series of ramps
% while undergoing intermittent collisions. Spatial Contact Force blocks
% are used to model the friction and normal forces between every pair of
% geometries that may potentially come into contact during the simulation
% (e.g., between one of the car's wheels and a railing). Each Spatial
% Contact Force block is able to generate brief high-impact contact forces
% to model collisions, as well as sustained contact forces to model rolling
% and sliding.
%
% Copyright 2019-2023 The MathWorks, Inc.

%%
open_system('SpatialContactForceBus');
%%
close_system('SpatialContactForceBus');
