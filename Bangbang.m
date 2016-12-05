%
% MATLAB BANG! BANG!
%
% 2-player aim game. Inspired by the 1990 Windows game "Bang! Bang!".
% See: http://playdosgamesonline.com/bang-bang.html
% This verision was developed for teaching purposes, as an example of game
% developed using Matlab that does not require an extensive programming 
% knowledge. 
%
% This implementation uses:
% - Selective execution (if/switch-case statements).
% - Recursive execution (for/while loops).
% - Basic input/output (input/fprintf).
% - Functions.
% - Anonymous functions.
% - Random number generators.
% - Plotting tools (lines, patches, and markers).
%
% It does NOT use:
% - Struct variables.
% - Cell variables.
% - Objects
% - GUI commands.
%
% Developed by: Joan Aguilar Mayans
% E-mail:       joana1@uci.edu
% Date:         16-Nov-2016

%
% Main script -------------------------------------------------------------
%

function BANGBANG

% Clear
clc;
close all;

% Intro
disp('%%%%%%%%%%%%%%%%%%%%%%');
disp('% MATLAB BANG! BANG! %');
disp('%%%%%%%%%%%%%%%%%%%%%%');
fprintf('\n');

% Ask the user for a terrain
fprintf('Select scenery. Options are:\n');
fprintf(' (1) Flat terrain\n');
fprintf(' (2) Tilted terrain\n');
fprintf(' (3) A hill in the middle\n');
fprintf(' (4) A tall hill in the middle\n');
fprintf(' (5) Player A on top of a hill\n');
fprintf(' (6) Player B on top of a hill\n');
fprintf(' (7) Random choice\n');
fprintf('\n');
ground_sel = input('Your choice: ');
fprintf('\n');
if ground_sel == 7
    ground_sel = randi(6);
end

% Assign ground level as an anonymous function
rand1 = rand();
rand2 = rand();
switch ground_sel
    case 1
        ground = @(x) zeros(size(x));  % Flat ground
    case 2
        ground = @(x) (0.4*rand1 - 0.2)*x;  % Tilted terrain
    case 3
        ground = @(x) (5*rand1 + 10)*cos(2*pi*(x - (40*rand2 + 30))/100)  + (0.4*rand2 - 0.2)*x;  % A hill in the middle
    case 4
        ground = @(x) (10*rand1 + 20)*cos(2*pi*(x - (40*rand2 + 30))/100)  + (0.4*rand2 - 0.2)*x;  % A tall hill in the middle
    case 5
        ground = @(x) (30*rand1 + 20)./(1 + exp(0.2*(x - 50)));  % A on top of a hill
    case 6
        ground = @(x) (30*rand1 + 20)./(1 + exp(-0.2*(x - 50)));  % B on top of a hill
    otherwise
        ground = @(x) zeros(size(x));  % Flat ground
end

% Assign players position
pos_A = [0; ground(0)];
pos_B = [100; ground(100)];

% Ask the user for wind
fprintf('Select wind. Options are:\n');
fprintf(' (1) No wind\n');
fprintf(' (2) Gentle breeze\n');
fprintf(' (3) Strong breeze\n');
fprintf(' (4) Strong gusty breeze\n');
fprintf(' (5) Strong gale\n');
fprintf(' (6) Strong gusty gale\n');
fprintf(' (7) Hurricane mode\n');
fprintf(' (8) Random choice\n');
fprintf('\n');
wind_sel = input('Your choice: ');
fprintf('\n');
if wind_sel == 8
    wind_sel = randi(7);
end

% Assign wind level as an anonymous function
rand3 = rand();
rand4 = rand();
switch wind_sel
    case 1  % No wind
        nominal_wind = 0;
        gusts = 0;
    case 2  % Gentle breeze
        nominal_wind = sign(rand3 - 0.5)*(1.7*rand4 + 3.6);
        gusts = 0;
    case 3  % Strong breeze
        nominal_wind = sign(rand3 - 0.5)*(2.5*rand4 + 11.4);
        gusts = 0;
    case 4  % Strong gusty breeze
        nominal_wind = sign(rand3 - 0.5)*(2.5*rand4 + 11.4);
        gusts = 1;
    case 5  % Strong gale
        nominal_wind = sign(rand3 - 0.5)*(3.1*rand4 + 21.1);
        gusts = 0;
    case 6  % Strong gusty gale
        nominal_wind = sign(rand3 - 0.5)*(3.1*rand4 + 21.1);
        gusts = 1;
    case 7  % Hurricane
        nominal_wind = sign(rand3 - 0.5)*(3.9*rand4 + 33.1);
        gusts = 1;
    otherwise
        nominal_wind = 0;
        gusts = 0;
end
if ~gusts
    wind = @(t, rand_val) [nominal_wind; 0];
else
    wind = @(t, rand_val) [1.5*nominal_wind + 0.5*nominal_wind*sin(2*pi*(t/20 + rand_val)); 0];
end

% Tell the user the wind situation
fprintf('Wind report:\n');
if nominal_wind == 0
    fprintf('There is no wind.\n');
else
    if nominal_wind > 0
        direction = 'right';
    else
        direction = 'left';
    end
    if ~gusts
        fprintf('Constant wind of %.2f m/s towards the %s.\n', nominal_wind, direction);
    else
        fprintf('Gusty wind of %.2f m/s towards the %s.\n', nominal_wind, direction);
    end
end
fprintf('\n');

% Initiliaze turn to player A and winner flag to false
winner_flag = 0;
turn = 'A';

% Plot scenery
fig = plotscenery(pos_A, pos_B, ground);

% Keep playing turns until there is a winner
while ~winner_flag
    
    % Display turn information and ask the player for angle and velocity
    fprintf('It''s %c''s turn.\n', turn);
    angle = input('Select angle (degrees): ');  % In degrees
    vel = input('Select velocity (m/s):  ');  % In m/s
    
    % Simulate each shot
    switch turn
        case 'A'
            [winner_flag, trajectory, dt] = fire(pos_A, angle, vel, pos_B, ground, wind);
        case 'B'
            angle = 180 - angle;
            [winner_flag, trajectory, dt] = fire(pos_B, angle, vel, pos_A, ground, wind);
    end
    
    % Plot the shot
    plottrajectory(fig, trajectory, dt, turn);
    
    % If there isn't a winner, tell the player that he/she missed and
    % change the turn variable
    if ~winner_flag
        fprintf('You missed!\n\n');
        switch turn
            case 'A'
                turn = 'B';
            case 'B'
                turn = 'A';
        end
    end
end

% After there is a winner, display it on screen.
fprintf('Player %c wins!\n', turn);

end

%
% Subfunctions ------------------------------------------------------------
%

function fig = plotscenery(pos_A, pos_B, ground)
% PLOTSCENERY generates a figure with the basic landscape for the game
% Bang! Bang! It places player A at pos_A, player B at pos_B and draws the
% ground using function ground.
% fig: figure handle where the scenery is plotted.
% pos_A: vector containing the x and y coordinates for player A.
% pos_B: vector containing the x and y coordinates for player B.
% ground: a function handle that defines the ground level (if y >
% ground(x), y is above the ground at coordinate x.
%
% Created by: Joan Aguilar Mayans
% E-mail:     joana1@uci.edu
% Date:       17-Nov-2016

% Create figure
fig = figure;
hold on;

% Generate axis
delta_hor = abs(pos_B(1) - pos_A(1));  % Horizontal distance between players
delta_ver = abs(pos_B(2) - pos_A(2));  % Vertical distance between players
max_delta = max(delta_hor, delta_ver);
mid = (pos_A + pos_B)/2;  % Middle point between the two players.
min_x_axis = mid(1) - max_delta/2 - 10;
max_x_axis = mid(1) + max_delta/2 + 10;
min_y_axis = mid(2) - max_delta/2 - 10;
max_y_axis = mid(2) + max_delta/2 + 10;
axis([min_x_axis, max_x_axis, min_y_axis, max_y_axis]);
axis square;

% Generate horizon points
horiz_x = linspace(min_x_axis, max_x_axis, 100);
horiz_y = ground(horiz_x);

% Plot ground and sky
patch([horiz_x, max_x_axis, min_x_axis], [horiz_y, min_y_axis, min_y_axis], [0, 0.6, 0]);
patch([horiz_x, max_x_axis, min_x_axis], [horiz_y, max_y_axis, max_y_axis], [0, 1, 1]);

% Plot players as markers
plot(pos_A(1), pos_A(2), 'Marker','s', 'MarkerFaceColor', 'red', 'MarkerSize', 20);
text(pos_A(1) - 5, pos_A(2) + 5, 'A');
plot(pos_B(1), pos_B(2), 'Marker','s', 'MarkerFaceColor', 'blue', 'MarkerSize', 20);
text(pos_B(1) + 5, pos_B(2) + 5, 'B');

end

function [hit, traj, dt] = fire(your_pos, angle, lin_vel, enemy_pos, ground, wind)
% FIRE computes the trajectory of a shot and checks if it hits the enemy.
% hit: either true or false depending if the shot ends up close enough to 
% enemy_pos.
% traj: an array containing points equally spaced in time from the shot
% your_pos: a column vector containing the position of the player that is 
% shooting.
% dt: time step used to compute the trajectory.
% angle: angle of the shot.
% vel: velocity of the shot.
% enemy_pos: a column vector containing the position of the player that is 
% being shot.
% ground: a function that defines the ground level (m) as a function of x 
% (if y > ground(x), y is above the ground at coordinate x).
% wind: a function that defines the wind speed (m/s) as a function of time 
% and a random parameter between 0 and 1
%
% Created by: Joan Aguilar Mayans
% E-mail:     joana1@uci.edu
% Date:       17-Nov-2016

% Note: variable allocation can be made more efficient.

% Define simulation parameters
m = 10;  % mass in kg;
r = 0.1;  % radius in m
rho_air = 1.225;  % air density at sea level in kg/m3
cd = 0.47;  % Sphere drag coefficient
A = pi*r^2;  % Reference area (m2)
g = [0; -9.81];  % Gravity vector in m/s2

% Define integration time step
dt = 0.01;

% Define wind random parameter
rand_par = rand();

% Initialize
t = 0;
traj = your_pos;
vel_vec = [lin_vel*cosd(angle);
           lin_vel*sind(angle)];

% Compute trajectory
while (traj(2,end) >= ground(traj(1,end)))
    
    % Compute current acceleration, velocity, and position
    u = wind(t, rand_par) - vel_vec(:,end);  % Flow velocity
    F_d = 0.5*rho_air*cd*A*norm(u)*u;  % Aerodynamic drag force
    acc = (m*g + F_d)/m;  % Compute acceleration
    vel = vel_vec(:,end);  % Compute velocity
    pos = traj(:,end); % Compute position
    
    % Compute next step velocity and position (semi-implicit Euler)
    vel_plus = vel + acc*dt;
    pos_plus = pos + vel_plus *dt;
    
    % Update time
    t = t + dt;
    
    % Store
    vel_vec = [vel_vec, vel_plus];
    traj = [traj, pos_plus];
end

% Check if there is a hit
hit = norm(traj(:,end) - enemy_pos) < 1;
end

function plottrajectory(fig, traj, dt, turn)
% PLOTTRAJECTORY plots the trajectory defined by traj on the scenery. The
% color depennds on the variable turn.
% fig: figure handle where the scenery is plotted.
% traj: an array containing the points of the trajectory to plot.
% dt: time step used to compute the trajectory.
% turn: is either A or B and determines the color in which the trajectory
% is plotted (A is red and B is blue).
%
% Created by: Joan Aguilar Mayans
% E-mail:     joana1@uci.edu
% Date:       17-Nov-2016

% Get number of points in the trajectory
N = size(traj, 2);

% Get color
if turn == 'A'
    color = [1, 0, 0];
elseif turn == 'B'
    color = [0, 0, 1];
end

% Plot
figure(fig);
for n = 2:N
    tic;
    plot([traj(1, n - 1), traj(1, n)], [traj(2, n - 1), traj(2, n)], '-', 'Color', color);
    pause(dt - toc);
end

% Emphasize last point
plot(traj(1, end), traj(2, end), 'Marker', 'o', 'MarkerFaceColor', 'black', 'MarkerSize', 20);

end
