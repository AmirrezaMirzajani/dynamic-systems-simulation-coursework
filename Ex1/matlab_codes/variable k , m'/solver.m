% main_simulation_timevar.m
% Full simulation with time-varying k1(t) and m'(t), including parameter plots

clear; clc; close all;

% Load parameters
p = params_timevar();

% Time settings
t_end = 10;
dt = 0.005;
t_span = 0:dt:t_end;
y0 = [0; 0; 0; 0];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Solve ODE
[T, Y] = ode45(@(t,y) equations_of_motion_timevar(t, y, p), t_span, y0, options);

% Extract results
x = Y(:,1);
x_dot = Y(:,2);
theta = Y(:,3);
theta_dot = Y(:,4);

% Compute accelerations
x_ddot = zeros(size(T));
theta_ddot = zeros(size(T));
for i = 1:length(T)
    [~, x_ddot(i), theta_ddot(i)] = equations_of_motion_with_accel_timevar(T(i), Y(i,:)', p);
end

% Compute time-varying parameters for plotting
m_prime_vals = zeros(size(T));
k1_vals = zeros(size(T));
for i = 1:length(T)
    m_prime_vals(i) = p.m_prime(T(i));
    k1_vals(i) = p.k1(T(i));
end

% ----- Plotting (black color, Times New Roman, no titles) -----
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
figWidth = 560;
figHeight = 420;
myColor = [0 0 0];

% 1. x(t)
figure('Name', 'x(t) time-varying', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T, x, 'Color', myColor, 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$x(t)$ [m]', 'Interpreter', 'latex', 'FontSize', 10);
grid on; box off; set(gca, 'FontSize', 10);

% 2. x_dot(t)
figure('Name', 'xdot(t) time-varying', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T, x_dot, 'Color', myColor, 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$\dot{x}(t)$ [m/s]', 'Interpreter', 'latex', 'FontSize', 10);
grid on; box off; set(gca, 'FontSize', 10);

% 3. x_ddot(t)
figure('Name', 'xddot(t) time-varying', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T, x_ddot, 'Color', myColor, 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$\ddot{x}(t)$ [m/s$^2$]', 'Interpreter', 'latex', 'FontSize', 10);
grid on; box off; set(gca, 'FontSize', 10);

% 4. theta(t)
figure('Name', 'theta(t) time-varying', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T, theta, 'Color', myColor, 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$\theta(t)$ [rad]', 'Interpreter', 'latex', 'FontSize', 10);
grid on; box off; set(gca, 'FontSize', 10);

% 5. theta_dot(t)
figure('Name', 'thetadot(t) time-varying', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T, theta_dot, 'Color', myColor, 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$\dot{\theta}(t)$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 10);
grid on; box off; set(gca, 'FontSize', 10);

% 6. theta_ddot(t)
figure('Name', 'thetaddot(t) time-varying', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T, theta_ddot, 'Color', myColor, 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$\ddot{\theta}(t)$ [rad/s$^2$]', 'Interpreter', 'latex', 'FontSize', 10);
grid on; box off; set(gca, 'FontSize', 10);

% 7. Phase portrait (x vs theta)
figure('Name', 'Phase portrait time-varying', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(x, theta, 'Color', myColor, 'LineWidth', 1);
xlabel('$x$ [m]', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('$\theta$ [rad]', 'Interpreter', 'latex', 'FontSize', 10);
grid on; box off; set(gca, 'FontSize', 10);

% 8. m'(t) parameter variation
figure('Name', 'm_prime(t)', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T, m_prime_vals, 'Color', myColor, 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$m''(t)$ [kg]', 'Interpreter', 'latex', 'FontSize', 10);
grid on; box off; set(gca, 'FontSize', 10);

% 9. k1(t) parameter variation
figure('Name', 'k1(t)', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T, k1_vals, 'Color', myColor, 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$k_1(t)$ [N/m]', 'Interpreter', 'latex', 'FontSize', 10);
grid on; box off; set(gca, 'FontSize', 10);

disp('Time-varying simulation completed. Results saved in workspace.');