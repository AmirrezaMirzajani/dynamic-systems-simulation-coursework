% run_simulation_case3_combined.m
% Forced harmonic excitation: F(t)=F0*sin(omega*t), T(t)=T0*sin(omega*t)
% Three frequencies in one figure for x(t) and one figure for theta(t)

clear; close all; clc;

%% Load parameters
p = parameters();

%% Excitation parameters
F0 = 10;                    % N
T0 = 5;                     % N.m
omega_list = [0.5, 1.2, 2.0];   % rad/s (three different frequencies)

%% Time span (long enough to see steady state)
tspan = [0, 30];

%% Define colors for different frequencies
colors = {'b', 'r', 'g'};
lineStyles = {'-', '-', '-'};

%% Storage for results
results = cell(length(omega_list), 1);

%% Solve ODE for each frequency
for i = 1:length(omega_list)
    w = omega_list(i);
    
    % Force functions for this omega
    F_fun = @(t) F0 * sin(w * t);
    T_fun = @(t) T0 * sin(w * t);
    
    % Zero initial conditions
    IC = [0; 0; 0; 0];
    
    % Solve ODE
    [t, y] = ode45(@(t,y) dynamics_nonlinear(t,y,p,F_fun,T_fun), tspan, IC);
    
    % Store results
    results{i}.t = t;
    results{i}.x = y(:,1);
    results{i}.theta = y(:,3);
    results{i}.omega = w;
end

%% Figure 1: Displacement x(t) for all three frequencies in one plot
figure;
hold on;
for i = 1:length(omega_list)
    plot(results{i}.t, results{i}.x, 'Color', colors{i}, 'LineWidth', 1.2, 'LineStyle', lineStyles{i});
end
hold off;
grid on;
xlabel('Time (s)');
ylabel('x(t) (m)');
legend({'\omega = 0.5 rad/s', '\omega = 1.2 rad/s', '\omega = 2.0 rad/s'}, 'Location', 'best');
% No title (you will add manually in Persian)

%% Figure 2: Angle theta(t) for all three frequencies in one plot
figure;
hold on;
for i = 1:length(omega_list)
    plot(results{i}.t, results{i}.theta, 'Color', colors{i}, 'LineWidth', 1.2, 'LineStyle', lineStyles{i});
end
hold off;
grid on;
xlabel('Time (s)');
ylabel('\theta(t) (rad)');
legend({'\omega = 0.5 rad/s', '\omega = 1.2 rad/s', '\omega = 2.0 rad/s'}, 'Location', 'best');
% No title

%% Optional: Display steady-state amplitude comparison
fprintf('Steady-state amplitude comparison (last 10 seconds):\n');
for i = 1:length(omega_list)
    t = results{i}.t;
    x = results{i}.x;
    theta = results{i}.theta;
    idx = t >= (max(t)-10);
    amp_x = (max(x(idx)) - min(x(idx))) / 2;
    amp_theta = (max(theta(idx)) - min(theta(idx))) / 2;
    fprintf('ω = %.1f rad/s: amp_x = %.4f m, amp_theta = %.4f rad\n', ...
        omega_list(i), amp_x, amp_theta);
end