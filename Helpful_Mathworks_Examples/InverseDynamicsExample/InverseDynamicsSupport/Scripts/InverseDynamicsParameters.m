% Color and density values
color = simmechanics.demohelpers.colors;
density = simmechanics.demohelpers.densities;

% Time constant
motion_tc = 0.02;

% Robot and Workpiece parameters
robot.base.height = 100; %cm
robot.base.radius = 12; %cm

robot.torso.height = 90; %cm
robot.torso.radius = 9; %cm

robot.head.length = 30; %cm
robot.head.width = 24; %cm
robot.head.height = 20; %cm

robot.lowerarmpin.radius = 1; %cm
robot.lowerarmpin.height = 3; %cm

robot.shoulder.offset = 15; %cm

robot.upperarm.length = 70; %cm
robot.upperarm.iradius = 8; %cm
robot.upperarm.oradius = 6; %cm
robot.upperarm.thickness = 6; %cm

robot.middlearm.length = 30; %cm
robot.middlearm.radius = 6; %cm
robot.middlearm.width = 10; %cm
robot.middlearm.thickness = 5; %cm

robot.lowerarm.length = 50; %cm
robot.lowerarm.radius = 3; %cm
robot.lowerarm.width = 10; %cm
robot.lowerarm.thickness = 5; %cm

robot.ujoint.plate.length = 55; %mm
robot.ujoint.plate.width = 40; %mm
robot.ujoint.plate.thickness = 10; %mm
robot.ujoint.axle.offset = 30; %mm
robot.ujoint.axle.radius = 10; %mm
robot.ujoint.crosspiece.length = 64; %mm
robot.ujoint.crosspiece.radius = 10; %mm

robot.gripper.basethickness = 10; %mm
robot.gripper.baseradius = 18; %mm
robot.gripper.tipradius = 10; %mm
robot.gripper.clawlength = 50; %mm
robot.gripper.tipdepth = 15; %mm
robot.gripper.channelradius = 5; %mm

robot.tool.length = 15; %cm
robot.tool.shaftradius = 0.5; %cm
robot.tool.tiplength = 1.5; %cm
robot.tool.tiporadius = 0.3; %cm
robot.tool.tipiradius = 1.2; %cm

robot.wrist.offset = 6; %cm

workpiece.distance = 120; %cm
workpiece.height = 150; %cm
workpiece.glyphsep = 40; %cm
workpiece.plate.length = 140; %cm
workpiece.plate.width = 50; %cm
workpiece.plate.thickness = 2; %cm
workpiece.beadsize = 0.5; %cm
workpiece.plus.armlength = 20; %cm
workpiece.circle.radius = 10; %cm
workpiece.star.radius = 10; %cm

floor.length = 50; %cm
floor.width = 50; %cm
floor.thickness = 4; %cm

% Tool tip home location in robot reference frame.
%
% The 5.0 cm offset in the x-coordinate is necessary to avoid a singularity
% in the robot's kinematics.
robot.tool.location = [ ...
    robot.upperarm.length + robot.wrist.offset / 2 + ...
    robot.tool.length + 5.0; ...
    -robot.shoulder.offset; ...
    robot.base.height + robot.head.height / 2 + ...
    robot.middlearm.length + robot.lowerarm.length + ...
    robot.wrist.offset / 2 ...
    ]; %cm

% Workpiece location in robot reference frame
robot.workpiece.location = [workpiece.distance 0 workpiece.height]'; %cm
robot.workpiece.rotation = [0 0 -1; -1 0 0; 0 1 0];

% Convert from robot ref coordinates to workpiece coordinates
workpiece.tool.home = robot.workpiece.rotation' * ...
                      (robot.tool.location - robot.workpiece.location); %cm
workpiece.tool.zoffset = 6; %cm
