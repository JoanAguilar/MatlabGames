%
% MATLAB MOON LANDER
%
% 2D Moon landing simulation game. Inspired by the 1990 Windows game 
% "Lander". See for instance:
% https://play.google.com/store/apps/details?id=com.pilot51.lander
% The player needs to land the lunar module safely on the Moon using a
% limitted amount of fuel. This game was developed as an example of a 
% simple Matlab game and to test different Matlab capabilities.
%
% Developed by: Joan Aguilar Mayans
% E-mail:       joana1@uci.edu
% Date:         21-Nov-2016

%
% Main script -------------------------------------------------------------
%

function Lander

% Clear
clear;
clc;
close all;

% Check if we are running Octave or Matlab
is_octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

% Intro
disp('%%%%%%%%%%%%%%%%%');
disp('% MATLAB LANDER %');
disp('%%%%%%%%%%%%%%%%%');
fprintf('\n');
fprintf('Your mission is to land the Apollo Lunar lander safely on the moon.\n');
fprintf('\n');

% Ask the user for control options
fprintf('Select controls. Options are:\n');
fprintf(' (1) Arrow keys:\n');
fprintf('      - Left arrow: RCS thrusters counter-clockwise.\n');
fprintf('      - Up arrow: Main thruster.\n');
fprintf('      - Right arrow: RCS thrusters clockwise.\n');
fprintf(' (2) Numpad:\n');
fprintf('      - 4: RCS thrusters counter-clockwise.\n');
fprintf('      - 7: RCS thrusters counter-clockwise + main thruster.\n');
fprintf('      - 8: Main thruster.\n');
fprintf('      - 9: RCS thrusters clockwise + main thruster.\n');
fprintf('      - 6: RCS thrusters clockwise.\n');
fprintf(' (3) AWSD:\n');
fprintf('      - A: RCS thrusters counter-clockwise.\n');
fprintf('      - Q: RCS thrusters counter-clockwise + main thruster.\n');
fprintf('      - W: Main thruster.\n');
fprintf('      - E: RCS thrusters clockwise + main thruster.\n');
fprintf('      - D: RCS thrusters clockwise.\n');
fprintf('\n');
controls_sel = input('Your choice: ');
fprintf('\n');

% Ask the user for a difficulty
fprintf('Select difficulty. Options are:\n');
fprintf(' (1) Easy:\n');
fprintf('      - 60 seconds left of fuel.\n');
fprintf('      - 10 m/s of max vertical velocity.\n');
fprintf('      - 5 m/s of max horizontal velocity.\n');
fprintf('      - 30 degrees of max terrain angle.\n');
fprintf('      - 30 degrees of max lander-terrain angle.\n');
fprintf(' (2) Medium:\n');
fprintf('      - 40 seconds left of fuel.\n');
fprintf('      - 5 m/s of max vertical velocity.\n');
fprintf('      - 2 m/s of max horizontal velocity.\n');
fprintf('      - 10 degrees of max terrain angle.\n');
fprintf('      - 10 degrees of max lander-terrain angle.\n');
fprintf(' (3) Hard:\n');
fprintf('      - 20 seconds left of fuel.\n');
fprintf('      - 2 m/s of max vertical velocity.\n');
fprintf('      - 1 m/s of max horizontal velocity.\n');
fprintf('      - 5 degrees of max terrain angle.\n');
fprintf('      - 5 degrees of max lander-terrain angle.\n');
fprintf(' (4) Impossible:\n');
fprintf('      - 10 seconds left of fuel.\n');
fprintf('      - 1 m/s of max vertical velocity.\n');
fprintf('      - 0.5 m/s of max horizontal velocity.\n');
fprintf('      - 2 degrees of max terrain angle.\n');
fprintf('      - 2 degrees of max lander-terrain angle.\n');
fprintf('\n');
diff_sel = input('Your choice: ');
fprintf('\n');

% Assign difficulty level variables
switch diff_sel
    case 1
        fuel_time = 60;
        max_vel_y = 10;
        max_vel_x = 5;
        max_terrain_angle = 30*pi/180;
        max_lander_angle = 30*pi/180;
    case 2
        fuel_time = 40;
        max_vel_y = 5;
        max_vel_x = 2;
        max_terrain_angle = 10*pi/180;
        max_lander_angle = 10*pi/180;
    case 3
        fuel_time = 20;
        max_vel_y = 2;
        max_vel_x = 1;
        max_terrain_angle = 5*pi/180;
        max_lander_angle = 5*pi/180;
    case 4
        fuel_time = 10;
        max_vel_y = 1;
        max_vel_x = 0.5;
        max_terrain_angle = 2*pi/180;
        max_lander_angle = 2*pi/180;
    otherwise
        fuel_time = 60;
        max_vel_y = 10;
        max_vel_x = 5;
        max_terrain_angle = 30*pi/180;
        max_lander_angle = 30*pi/180;
end

% Landing site dimensions
land_dimensions = [-50 50 0 100];

% Generate stars positions
N_stars = 200;
stars = [(land_dimensions(2) - land_dimensions(1))*rand(1,N_stars) + land_dimensions(1);
         (land_dimensions(4) - land_dimensions(3))*rand(1,N_stars) + land_dimensions(3)];

% Definie lander initial state
pos = [30; 100];  % Position (m)
R = diag([1 1]); % Rotation matrix
vel = [0; 0];    % Velocity (m/s)
w = 0;           % Angular velocity (rad/s)
T_B = [0; 0];    % Thrust in body frame (N)
M_G = 0;         % Applied torque (at the center of mass, Nm)
m_dry = 4284;    % Dry mass (kg)
m_fuel = 0.5*10463;  % Fuel mass (kg)
m_RCS = 0.75*287;     % RCS fuel mass (kg)
mdot_fuel = -m_fuel/fuel_time;  % Fuel flow (kg/s)
mdot_RCS = -m_RCS/(2*fuel_time);  % RCS fuel flow (kg/s)
contact_point_B = [0; -4.2];  % Contact point with the ground in body coordinates (m)

% Define lunar gravity
g = [0; -1.62];

% Define figure
fig = figure;
hold on;

% Plot scene
plotscene(fig, land_dimensions, stars, @ground, pos, R, T_B, M_G);

% Display information on the command window
dispinfo(@ground, pos, R, vel, w, contact_point_B, m_dry, m_fuel, m_RCS, is_octave);

% Initialize values for simulation
% key = 'downarrow';
landed = 0;  % 0 if the lander is flying
             % 1 if it has landed 
             % 2 if it has crashed
             % 3 if the user quits
dt = 0.1;  % Time step
while landed == 0
    
    % Set counter
    tic;
    
    % Get the key that is being pressed
    set(fig, 'KeyPressFcn', @keypresscallback, 'KeyReleaseFcn', @keyreleasecallback);
    key = getappdata(0, 'key');
    if isempty(key)
        key = 'downarrow';
    end
    
    % Exit if the user hit escape
    if strcmp(key, 'escape')
        landed = 3;
        break;
    end
        

    % Compute thrust and recompute mass
    if m_fuel > 0
        switch controls_sel
            case 2
                switch key
                    case {'numpad7', 'numpad8', 'numpad9'}
                        T_B = [0; 45040];
                        m_fuel = m_fuel + mdot_fuel*dt;
                    otherwise
                        T_B = [0; 0];
                end
            case 3
                switch key
                    case {'q', 'w', 'e'}
                        T_B = [0; 45040];
                        m_fuel = m_fuel + mdot_fuel*dt;
                    otherwise
                        T_B = [0; 0];
                end     
            otherwise
                switch key
                    case 'uparrow'
                        T_B = [0; 45040];
                        m_fuel = m_fuel + mdot_fuel*dt;
                    otherwise
                        T_B = [0; 0];
                end
        end
    else
        T_B = [0; 0];
    end
    
    % Compute torque and recompute mass
    if m_RCS > 0
        switch controls_sel
            case 2
                switch key
                    case {'numpad4', 'numpad7'}
                        M_G = 2*440*4.29;
                        m_RCS = m_RCS + mdot_RCS*dt;
                    case {'numpad9', 'numpad6'}
                        M_G = -2*440*4.29;
                        m_RCS = m_RCS + mdot_RCS*dt;
                    otherwise
                        M_G = 0;
                end
            case 3
                switch key
                    case {'a', 'q'}
                        M_G = 2*440*4.29;
                        m_RCS = m_RCS + mdot_RCS*dt;
                    case {'e', 'd'}
                        M_G = -2*440*4.29;
                        m_RCS = m_RCS + mdot_RCS*dt;
                    otherwise
                        M_G = 0;
                end
            otherwise
                switch key
                    case 'leftarrow'
                        M_G = 2*440*4.29;
                        m_RCS = m_RCS + mdot_RCS*dt;
                    case 'rightarrow'
                        M_G = -2*440*4.29;
                        m_RCS = m_RCS + mdot_RCS*dt;
                    otherwise
                        M_G = 0;
                end
        end
    else
        M_G = 0;
    end
    
    % Recompute mass
    if m_fuel < 0
        m_fuel = 0;
    end
    if m_RCS < 0
        m_RCS = 0;
    end
    m_total = m_dry + m_fuel + m_RCS;
    
    % Recompute moment of inertia (as a sphere of radius = 3 m)
    I_G = 2*m_total*3^2/5;
    
    % Integrate velocities
    w = w + (1/I_G)*M_G*dt;
    vel = vel + (g + (1/m_total)*(R*T_B))*dt;
    
    % Integrate position and orientation
    th = w*dt;
    Rstep = [cos(th), -sin(th);
             sin(th),  cos(th)];
    R = Rstep*R;
    pos = pos + vel*dt;
    
    % Replot scene
    % Plot scene
    plotscene(fig, land_dimensions, stars, @ground, pos, R, T_B, M_G);

    % Update information on the command window
    dispinfo(@ground, pos, R, vel, w, contact_point_B, m_dry, m_fuel, m_RCS, is_octave);
    
    % Check if flying, landed, or crashed
    contact_point = R*contact_point_B + pos;
    ground_dist = contact_point(2) - ground(contact_point(1));
    if ground_dist > 0
        landed = 0;
    else
        [~, dydx] = ground(pos(1));
        terrain_angle = atan2(dydx, 1);
        lander_angle = sign(R(2,2))*acos(trace(R)/2);
        if -vel(2) < max_vel_y && abs(vel(1)) < max_vel_x && abs(terrain_angle) < max_terrain_angle && abs(lander_angle - terrain_angle) < max_lander_angle
            landed = 1;
        else
            landed = 2;
        end
    end
       
    % Pause if needed
    pause(dt - toc);
end

fprintf('\n');
if landed == 1
    fprintf('You landed!\n');
elseif landed == 2
    fprintf('You crashed!\n');
elseif landed == 3
    fprintf('Bye!\n');
end

end

% % Flat ground
% function [y, dydx] = ground(x)
% % GROUND returns the ground level (y) and the gradient (dydx). y values
% % should be at least 1.
% 
% y = ones(size(x));
% dydx = zeros(size(x));
% 
% end

% Crater
function [y, dydx] = ground(x)
% GROUND returns the ground level (y) and the gradient (dydx). y values
% should be at least 1.

a = 10;
b = pi/40;
c = 1;

y = a.^(-cos(b*x)) + c;
dydx = (b*sin(b*x).*log(a))./a.^cos(b*x);

end

function plotscene(fig, full_dimensions, stars, ground, pos, R, T, M)
% PLOTSCENE plots a region of the scene of the lunar landing. The region
% can be of 100x100 (full scene), 50x50, 20x20, or 10x10.
% fig: a value containing the figure number.
% dimensions: a vector containing the dimensions of the landing site using
% format [xmin xmax ymin ymax].
% stars: 2xN array containing the position of stars.
% ground: a function that returns the ground level at a specific point and
% the gradient at that point.
% pos: a vector containing the current position of the center of mass of
% the lander (m).
% R: a rotation matrix containing the orientation of the lander.
% T: a vector containing the thrust in the body frame (N).
% M: a scalar containing the torque applied about the center of mass (Nm).

% Select which region to plot
if pos(2) > 50
    region = full_dimensions;
else
    if pos(2) > 20
        scene_side = 50;
    elseif pos(2) > 10
        scene_side = 20;
    else
        scene_side = 10;
    end
    region = [scene_side*floor(pos(1)/scene_side), scene_side*floor(pos(1)/scene_side) + scene_side, 0, scene_side];
end

% Figure out which starts to plot
plotted_stars = (stars(1,:) < region(2) & stars(1,:) > region(1) & stars(2,:) < region(4) & stars(2,:) > region(3));

% Compute ground level
x_ground = linspace(region(1), region(2), 100);
y_ground = ground(x_ground);

% Select figure
figure(fig);
cla;

% Plot sky
patch([region(1), region(2), region(2), region(1)], [region(3), region(3), region(4), region(4)], 'black');

% Plot stars
plot(stars(1, plotted_stars), stars(2, plotted_stars), 'w.');

% Plot lander
plotlander(fig, pos, R, T, M);

% Plot ground
patch([x_ground, region(2), region(1)], [y_ground, region(3), region(3)], [.6 .6 .6]);

% Set axes
axis(region);
axis square;

% Force draw
drawnow;

end

function plotlander(fig, pos, R, T, M)
% PLOTLANDER plots a representation of the Lunar Lander on the scene. It
% may also plot the flames coming from the different thrusters.

% Vertices of the ascent stage in body frame
ah = 2.832;  % Ascent stage height (m)
aw = 4.29;   % Ascent stage width (m)
ascent_vertices_B = [-aw/2, aw/2, aw/2,         0,     -aw/2;
                     -1,   -1,    (2/3)*ah - 1, ah - 1, (2/3)*ah - 1];
N_ascent_vertices = size(ascent_vertices_B, 2);

% Vertices of th descent stage in body frame
dh = 3.231;  % Descent stage height (m)
dw = 9.4;   % Ascent stage width (m)
descent_vertices_B = [-dw/2,    0,         dw/2,    aw/2, -aw/2;
                      -dh - 1, -dh/2 - 1, -dh - 1, -1,    -1];
N_descent_vertices = size(descent_vertices_B, 2);

% Plot thruster if needed
if T(2) > 0
    thrust_vertices_B = (T(2)/45000)*[-1,  0,  1, 0;  % 45000 is the normalizing/max thrust
                                      -2, -4, -2, 0];
	thrust_vertices_B(2,:) = thrust_vertices_B(2,:) - dh;
	N_thrust_vertices = size(thrust_vertices_B, 2);
end

% Plot RCS thrusters if needed
if M ~= 0
    RCS_vertices_B = (M/45000)*[-1,  0,  1, 0;  % 45000 is the normalizing/max thrust
                                -2, -4, -2, 0];
    RCS_vertices_BR = RCS_vertices_B;
    RCS_vertices_BL = -RCS_vertices_B;
    RCS_vertices_BR(1,:) = RCS_vertices_BR(1,:) + aw/2 + 0.1;
    RCS_vertices_BL(1,:) = RCS_vertices_BL(1,:) - aw/2 - 0.1;
    N_RCS_vertices = size(RCS_vertices_B, 2);
end

% Rotate vertices
ascent_vertices = zeros(size(ascent_vertices_B));
for n = 1:N_ascent_vertices
    ascent_vertices(:,n) = R*ascent_vertices_B(:,n);
end
descent_vertices = zeros(size(descent_vertices_B));
for n = 1:N_descent_vertices
    descent_vertices(:,n) = R*descent_vertices_B(:,n);
end
if T(2) > 0
    thrust_vertices = zeros(size(thrust_vertices_B));
    for n = 1:N_thrust_vertices
        thrust_vertices(:,n) = R*thrust_vertices_B(:,n);
    end
end
if M ~= 0
    RCS_vertices_R = zeros(size(RCS_vertices_BR));
    RCS_vertices_L = zeros(size(RCS_vertices_BL));
    for n = 1:N_RCS_vertices
        RCS_vertices_R(:,n) = R*RCS_vertices_BR(:,n);
        RCS_vertices_L(:,n) = R*RCS_vertices_BL(:,n);
    end
end
    

% Translate vertices
ascent_vertices(1,:) =  ascent_vertices(1,:) +  pos(1);
ascent_vertices(2,:) =  ascent_vertices(2,:) +  pos(2);
descent_vertices(1,:) = descent_vertices(1,:) + pos(1);
descent_vertices(2,:) = descent_vertices(2,:) + pos(2);
if T(2) > 0
    thrust_vertices(1,:) = thrust_vertices(1,:) + pos(1);
    thrust_vertices(2,:) = thrust_vertices(2,:) + pos(2);
end
if M ~= 0
    RCS_vertices_R(1,:) = RCS_vertices_R(1,:) + pos(1);
    RCS_vertices_R(2,:) = RCS_vertices_R(2,:) + pos(2);
    RCS_vertices_L(1,:) = RCS_vertices_L(1,:) + pos(1);
    RCS_vertices_L(2,:) = RCS_vertices_L(2,:) + pos(2);
end

% Plot lander
figure(fig);
patch(ascent_vertices(1,:),  ascent_vertices(2,:),  [.6 .6 .6]);
patch(descent_vertices(1,:), descent_vertices(2,:), [.8 .6 .1]);

% Plot thruster (if needed)
if T(2) > 0
    patch(thrust_vertices(1,:), thrust_vertices(2,:), [1, .8, 0]);
end

% Plot RCS thrusters (if needed)
if M ~= 0
    patch(RCS_vertices_R(1,:), RCS_vertices_R(2,:), [1, .8, 0]);
    patch(RCS_vertices_L(1,:), RCS_vertices_L(2,:), [1, .8, 0]);
end

end

function dispinfo(ground, pos, R, vel, w, contact_point_B, m_dry, m_fuel, m_RCS, is_octave)
% DISPINFO displays on the command window information about the lunar
% lander. It displays: position, velocity, distance to the ground, total
% mass of the lander, and fuel left.

% Compute angle
theta = -sign(R(1,2))*(180/pi)*acos(trace(R)/2);

% Compute distance to the ground
contact_point = R*contact_point_B + pos;
ground_dist = contact_point(2) - ground(contact_point(1));

% Clear command window;
clc;

% Display information
fprintf('               x          y\n');
fprintf('Position: %6.2f m   %6.2f m\n', pos(1), pos(2));
fprintf('Velocity: %6.2f m/s %6.2f m/s\n', vel(1), vel(2));
fprintf('Angle to horizon:   %7.2f degrees\n', theta);
fprintf('Angular velocity:   %7.2f rad/s\n', w);
fprintf('Ground dist:        %7.2f m\n', ground_dist);
fprintf('Thruster fuel left: %7.2f kg\n', m_fuel);
fprintf('RCS fuel left:      %7.2f kg\n', m_RCS);
fprintf('Total mass:         %7.2f kg\n', m_dry + m_fuel + m_RCS);
if is_octave
    fflush(stdout);
end

end



function keypresscallback(fig, evt)

setappdata(0, 'key', evt.Key);

end

function keyreleasecallback(fig, evt)

setappdata(0, 'key', []);

end