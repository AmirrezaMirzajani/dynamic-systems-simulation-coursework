% run_case2_2.m
clear; close all; clc;

p = parameters();

% Force functions: F(t)=0, T(t) is a pulse
F_fun = @(t) 0;
T_fun = @(t) (t <= 0.05) * 20;   % 20 N·m for 0.05 s

% Zero initial conditions
IC = [0; 0; 0; 0];
tspan = [0, 10];

% Solve
[t, y] = ode45(@(t,y) dynamics_nonlinear(t,y,p,F_fun,T_fun), tspan, IC);

% Extract results
x = y(:,1);
theta = y(:,3);

% Plot x(t)
figure;
plot(t, x, 'b-', 'LineWidth', 1.5);
grid on; xlabel('Time (s)'); ylabel('x(t) (m)');
legend('Scenario 2-2: Impact on pendulum');

% Plot theta(t)
figure;
plot(t, theta, 'r-', 'LineWidth', 1.5);
grid on; xlabel('Time (s)'); ylabel('\theta(t) (rad)');
legend('Scenario 2-2: Impact on pendulum');