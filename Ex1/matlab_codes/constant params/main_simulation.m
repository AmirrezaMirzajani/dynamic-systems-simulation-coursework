% main_simulation.m
% Solve the 2-DOF vibration system using ode45

clear; clc; close all;

% Load parameters
p = params();

% Simulation time
t_start = 0;
t_end = 10;
dt = 0.005;
t_span = t_start:dt:t_end;

% Initial conditions
y0 = [0.001; 0; 0; 0];

% ODE options
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Solve
[T, Y] = ode45(@(t,y) equations_of_motion(t, y, p), t_span, y0, options);

% Extract results
x = Y(:,1);
x_dot = Y(:,2);
theta = Y(:,3);
theta_dot = Y(:,4);

% Compute accelerations (x_ddot and theta_ddot) using equations of motion
x_ddot = zeros(size(T));
theta_ddot = zeros(size(T));
for i = 1:length(T)
    y_i = Y(i,:)';
    [~, x_ddot_i, theta_ddot_i] = equations_of_motion_with_accel(T(i), y_i, p);
    x_ddot(i) = x_ddot_i;
    theta_ddot(i) = theta_ddot_i;
end

% Save results
save('results.mat', 'T', 'x', 'x_dot', 'theta', 'theta_dot');

% Plot results
plot_results;