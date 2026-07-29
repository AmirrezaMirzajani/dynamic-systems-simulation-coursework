% sensitivity_time_frf.m
% Sensitivity analysis for k1, k2, c1:
% For each parameter, three values -> plot time responses (x, theta, x_ddot, theta_ddot)
% and frequency response (amplitude of x vs omega)

clear; clc; close all;

% ----- Base parameters (constant) -----
p_base.M = 2.0;
p_base.I_beam = 0.5;
p_base.m_prime = 0.5;
p_base.L_prime = 0.2;
p_base.e = 0.05;
p_base.L = 0.8;
p_base.m = 0.1;
p_base.omega_fixed = 20;      % fixed excitation for time response
p_base.c1 = 10;               % default
p_base.d = 0.1;
p_base.k1 = 800;
p_base.k2 = 800;

% Compute I_m_prime (constant)
A_s = p_base.L_prime^2;
A_c = pi * p_base.e^2;
I_m_prime = p_base.m_prime * ( (1/6)*p_base.L_prime^4 - 0.5*pi*p_base.e^4 ) / (A_s - A_c);
p_base.I_m_prime = I_m_prime;
p_base.I_total = p_base.I_beam + I_m_prime + p_base.m_prime * p_base.d^2;

% Time simulation settings
t_end = 5;          % simulation time (seconds)
dt = 0.005;
t_span = 0:dt:t_end;
y0 = [0;0;0;0];
options = odeset('RelTol',1e-6,'AbsTol',1e-8);
steady_samples = round(5 / dt);   % for FRF (last 5 sec)

% ----- Frequency response settings -----
omega_range = 1:0.5:50;    % rad/s

% ----- Define parameter sets -----
k1_vals = [600, 800, 1000];
k2_vals = [600, 800, 1000];
c1_vals = [5, 20, 50];

colors = {'b', 'r', 'k'};
styles = {'-.', '-', '--'};

% ----- Sensitivity for k1 -----
figure('Name', 'Sensitivity to k1 - Time responses', 'Position', [100 100 800 600]);
for idx = 1:3
    p = p_base;
    p.k1 = k1_vals(idx);
    p.k2 = 800;
    p.c1 = 10;
    p.I_total = p.I_beam + p.I_m_prime + p.m_prime * p.d^2;
    [T, x, theta, x_ddot, theta_ddot] = run_time_simulation(p, t_span, y0, options);
    
    subplot(2,2,1); hold on;
    plot(T, x, 'Color', colors{idx}, 'LineStyle', styles{idx}, 'LineWidth', 1);
    subplot(2,2,2); hold on;
    plot(T, theta, 'Color', colors{idx}, 'LineStyle', styles{idx}, 'LineWidth', 1);
    subplot(2,2,3); hold on;
    plot(T, x_ddot, 'Color', colors{idx}, 'LineStyle', styles{idx}, 'LineWidth', 1);
    subplot(2,2,4); hold on;
    plot(T, theta_ddot, 'Color', colors{idx}, 'LineStyle', styles{idx}, 'LineWidth', 1);
end
subplot(2,2,1); xlabel('Time (s)'); ylabel('$x(t)$ [m]', 'Interpreter','latex'); grid on; box off; legend('k1=600','k1=800','k1=1000','Location','best');
subplot(2,2,2); xlabel('Time (s)'); ylabel('$\theta(t)$ [rad]', 'Interpreter','latex'); grid on; box off;
subplot(2,2,3); xlabel('Time (s)'); ylabel('$\ddot{x}(t)$ [m/s$^2$]', 'Interpreter','latex'); grid on; box off;
subplot(2,2,4); xlabel('Time (s)'); ylabel('$\ddot{\theta}(t)$ [rad/s$^2$]', 'Interpreter','latex'); grid on; box off;
sgtitle('');

% FRF for k1
figure('Name', 'FRF sensitivity to k1', 'Position', [100 100 560 420]);
hold on;
for idx = 1:3
    p = p_base;
    p.k1 = k1_vals(idx);
    p.k2 = 800;
    p.c1 = 10;
    p.I_total = p.I_beam + p.I_m_prime + p.m_prime * p.d^2;
    [omega_range, amp_x] = compute_frf(p, omega_range, y0, options, steady_samples);
    plot(omega_range, amp_x, 'Color', colors{idx}, 'LineStyle', styles{idx}, 'LineWidth', 1);
end
xlabel('\omega (rad/s)', 'FontSize', 10);
ylabel('$|x|$ (m)', 'Interpreter','latex', 'FontSize', 10);
legend('k1=600','k1=800','k1=1000','Location','best');
grid on; box off; title('');

% ----- Sensitivity for k2 -----
figure('Name', 'Sensitivity to k2 - Time responses', 'Position', [100 100 800 600]);
for idx = 1:3
    p = p_base;
    p.k1 = 800;
    p.k2 = k2_vals(idx);
    p.c1 = 10;
    [T, x, theta, x_ddot, theta_ddot] = run_time_simulation(p, t_span, y0, options);
    
    subplot(2,2,1); hold on;
    plot(T, x, 'Color', colors{idx}, 'LineStyle', styles{idx}, 'LineWidth', 1);
    subplot(2,2,2); hold on;
    plot(T, theta, 'Color', colors{idx}, 'LineStyle', styles{idx}, 'LineWidth', 1);
    subplot(2,2,3); hold on;
    plot(T, x_ddot, 'Color', colors{idx}, 'LineStyle', styles{idx}, 'LineWidth', 1);
    subplot(2,2,4); hold on;
    plot(T, theta_ddot, 'Color', colors{idx}, 'LineStyle', styles{idx}, 'LineWidth', 1);
end
subplot(2,2,1); xlabel('Time (s)'); ylabel('$x(t)$ [m]', 'Interpreter','latex'); grid on; box off; legend('k2=600','k2=800','k2=1000','Location','best');
subplot(2,2,2); xlabel('Time (s)'); ylabel('$\theta(t)$ [rad]', 'Interpreter','latex'); grid on; box off;
subplot(2,2,3); xlabel('Time (s)'); ylabel('$\ddot{x}(t)$ [m/s$^2$]', 'Interpreter','latex'); grid on; box off;
subplot(2,2,4); xlabel('Time (s)'); ylabel('$\ddot{\theta}(t)$ [rad/s$^2$]', 'Interpreter','latex'); grid on; box off;
sgtitle('');

% FRF for k2
figure('Name', 'FRF sensitivity to k2', 'Position', [100 100 560 420]);
hold on;
for idx = 1:3
    p = p_base;
    p.k1 = 800;
    p.k2 = k2_vals(idx);
    p.c1 = 10;
    [omega_range, amp_x] = compute_frf(p, omega_range, y0, options, steady_samples);
    plot(omega_range, amp_x, 'Color', colors{idx}, 'LineStyle', styles{idx}, 'LineWidth', 1);
end
xlabel('\omega (rad/s)', 'FontSize', 10);
ylabel('$|x|$ (m)', 'Interpreter','latex', 'FontSize', 10);
legend('k2=600','k2=800','k2=1000','Location','best');
grid on; box off; title('');

% ----- Sensitivity for c1 -----
figure('Name', 'Sensitivity to c1 - Time responses', 'Position', [100 100 800 600]);
for idx = 1:3
    p = p_base;
    p.k1 = 800;
    p.k2 = 800;
    p.c1 = c1_vals(idx);
    [T, x, theta, x_ddot, theta_ddot] = run_time_simulation(p, t_span, y0, options);
    
    subplot(2,2,1); hold on;
    plot(T, x, 'Color', colors{idx}, 'LineStyle', styles{idx}, 'LineWidth', 1);
    subplot(2,2,2); hold on;
    plot(T, theta, 'Color', colors{idx}, 'LineStyle', styles{idx}, 'LineWidth', 1);
    subplot(2,2,3); hold on;
    plot(T, x_ddot, 'Color', colors{idx}, 'LineStyle', styles{idx}, 'LineWidth', 1);
    subplot(2,2,4); hold on;
    plot(T, theta_ddot, 'Color', colors{idx}, 'LineStyle', styles{idx}, 'LineWidth', 1);
end
subplot(2,2,1); xlabel('Time (s)'); ylabel('$x(t)$ [m]', 'Interpreter','latex'); grid on; box off; legend('c1=5','c1=20','c1=50','Location','best');
subplot(2,2,2); xlabel('Time (s)'); ylabel('$\theta(t)$ [rad]', 'Interpreter','latex'); grid on; box off;
subplot(2,2,3); xlabel('Time (s)'); ylabel('$\ddot{x}(t)$ [m/s$^2$]', 'Interpreter','latex'); grid on; box off;
subplot(2,2,4); xlabel('Time (s)'); ylabel('$\ddot{\theta}(t)$ [rad/s$^2$]', 'Interpreter','latex'); grid on; box off;
sgtitle('');

% FRF for c1
figure('Name', 'FRF sensitivity to c1', 'Position', [100 100 560 420]);
hold on;
for idx = 1:3
    p = p_base;
    p.k1 = 800;
    p.k2 = 800;
    p.c1 = c1_vals(idx);
    [omega_range, amp_x] = compute_frf(p, omega_range, y0, options, steady_samples);
    plot(omega_range, amp_x, 'Color', colors{idx}, 'LineStyle', styles{idx}, 'LineWidth', 1);
end
xlabel('\omega (rad/s)', 'FontSize', 10);
ylabel('$|x|$ (m)', 'Interpreter','latex', 'FontSize', 10);
legend('c1=5','c1=20','c1=50','Location','best');
grid on; box off; title('');

disp('Sensitivity analysis (time and FRF) completed.');

% ----- Local helper functions -----
function [T, x, theta, x_ddot, theta_ddot] = run_time_simulation(p, t_span, y0, options)
    p.Fv = @(t) p.m * p.e * p.omega_fixed^2 * sin(p.omega_fixed * t);
    p.Fh = @(t) p.m * p.e * p.omega_fixed^2 * cos(p.omega_fixed * t);
    [T, Y] = ode45(@(t,y) equations_of_motion_constant(t,y,p), t_span, y0, options);
    x = Y(:,1);
    theta = Y(:,3);
    x_ddot = zeros(size(T));
    theta_ddot = zeros(size(T));
    for i = 1:length(T)
        [~, x_ddot(i), theta_ddot(i)] = equations_of_motion_with_accel_constant(T(i), Y(i,:)', p);
    end
end

function [omega_range, amp_x] = compute_frf(p, omega_range, y0, options, steady_samples)
    amp_x = zeros(size(omega_range));
    for i = 1:length(omega_range)
        p.omega = omega_range(i);
        p.Fv = @(t) p.m * p.e * p.omega^2 * sin(p.omega * t);
        p.Fh = @(t) p.m * p.e * p.omega^2 * cos(p.omega * t);
        [~, Y] = ode45(@(t,y) equations_of_motion_constant(t,y,p), [0 30], y0, options);
        if length(Y) > steady_samples
            x_steady = Y(end-steady_samples:end,1);
            amp_x(i) = (max(x_steady) - min(x_steady)) / 2;
        else
            amp_x(i) = NaN;
        end
    end
end