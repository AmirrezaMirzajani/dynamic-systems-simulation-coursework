% compare_symmetric_asymmetric.m
% Compare system response (including accelerations) for symmetric and asymmetric springs

clear; clc; close all;

% Define common parameters (based on previous params.m)
p.M = 2.0;
p.I_beam = 0.5;
p.m_prime = 0.5;
p.L_prime = 0.2;
p.e = 0.05;
p.d = 0.1;
p.L = 0.8;
p.c1 = 10;
p.m = 0.1;
p.omega = 20;  % fixed excitation frequency for time response

% Compute I_m_prime and I_total (same for both cases)
A_s = p.L_prime^2;
A_c = pi * p.e^2;
p.I_m_prime = p.m_prime * ( (1/6)*p.L_prime^4 - 0.5*pi*p.e^4 ) / (A_s - A_c);
p.I_total = p.I_beam + p.I_m_prime + p.m_prime * p.d^2;

% Time settings
t_end = 10;
dt = 0.005;
t_span = 0:dt:t_end;
y0 = [0;0;0;0];
options = odeset('RelTol',1e-6,'AbsTol',1e-8);

% ----- Case 1: Symmetric springs (k1 = k2 = 800) -----
p.k1 = @(t) 800;
p.k2 = @(t) 800;
p.Fv = @(t) p.m * p.e * p.omega^2 * sin(p.omega * t);
p.Fh = @(t) p.m * p.e * p.omega^2 * cos(p.omega * t);

[T_sym, Y_sym] = ode45(@(t,y) equations_of_motion(t,y,p), t_span, y0, options);
x_sym = Y_sym(:,1);
x_dot_sym = Y_sym(:,2);
theta_sym = Y_sym(:,3);
theta_dot_sym = Y_sym(:,4);

% Compute accelerations for symmetric case
x_ddot_sym = zeros(size(T_sym));
theta_ddot_sym = zeros(size(T_sym));
for i = 1:length(T_sym)
    y_i = Y_sym(i,:)';
    [~, x_ddot_i, theta_ddot_i] = equations_of_motion_with_accel(T_sym(i), y_i, p);
    x_ddot_sym(i) = x_ddot_i;
    theta_ddot_sym(i) = theta_ddot_i;
end

% ----- Case 2: Asymmetric springs (k1=800, k2=500) -----
p.k1 = @(t) 800;
p.k2 = @(t) 500;
% Forces remain same
[T_asym, Y_asym] = ode45(@(t,y) equations_of_motion(t,y,p), t_span, y0, options);
x_asym = Y_asym(:,1);
x_dot_asym = Y_asym(:,2);
theta_asym = Y_asym(:,3);
theta_dot_asym = Y_asym(:,4);

% Compute accelerations for asymmetric case
x_ddot_asym = zeros(size(T_asym));
theta_ddot_asym = zeros(size(T_asym));
for i = 1:length(T_asym)
    y_i = Y_asym(i,:)';
    [~, x_ddot_i, theta_ddot_i] = equations_of_motion_with_accel(T_asym(i), y_i, p);
    x_ddot_asym(i) = x_ddot_i;
    theta_ddot_asym(i) = theta_ddot_i;
end

% ----- Plotting (symmetric: black dashed, asymmetric: red solid, linewidth 1) -----
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
figWidth = 560;
figHeight = 420;

% 1. Comparison of x(t)
figure('Name', 'Comparison x(t)', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T_sym, x_sym, 'k--', 'LineWidth', 1); hold on;
plot(T_asym, x_asym, 'r-', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$x(t)$ [m]', 'Interpreter', 'latex', 'FontSize', 10);
legend('Symmetric (k1=k2=800)', 'Asymmetric (k1=800,k2=500)', 'Location', 'best');
grid on; box off; set(gca, 'FontSize', 10);

% 2. Comparison of theta(t)
figure('Name', 'Comparison theta(t)', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T_sym, theta_sym, 'k--', 'LineWidth', 1); hold on;
plot(T_asym, theta_asym, 'r-', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$\theta(t)$ [rad]', 'Interpreter', 'latex', 'FontSize', 10);
legend('Symmetric', 'Asymmetric', 'Location', 'best');
grid on; box off; set(gca, 'FontSize', 10);

% 3. Comparison of x_ddot(t)
figure('Name', 'Comparison x_ddot(t)', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T_sym, x_ddot_sym, 'k--', 'LineWidth', 1); hold on;
plot(T_asym, x_ddot_asym, 'r-', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$\ddot{x}(t)$ [m/s$^2$]', 'Interpreter', 'latex', 'FontSize', 10);
legend('Symmetric', 'Asymmetric', 'Location', 'best');
grid on; box off; set(gca, 'FontSize', 10);

% 4. Comparison of theta_ddot(t)
figure('Name', 'Comparison theta_ddot(t)', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T_sym, theta_ddot_sym, 'k--', 'LineWidth', 1); hold on;
plot(T_asym, theta_ddot_asym, 'r-', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$\ddot{\theta}(t)$ [rad/s$^2$]', 'Interpreter', 'latex', 'FontSize', 10);
legend('Symmetric', 'Asymmetric', 'Location', 'best');
grid on; box off; set(gca, 'FontSize', 10);

% 5. Frequency response comparison
omega_range = 1:0.5:50;
amp_x_sym = zeros(size(omega_range));
amp_x_asym = zeros(size(omega_range));

% Symmetric FRF
p.k1 = @(t) 800; p.k2 = @(t) 800;
for i = 1:length(omega_range)
    p.omega = omega_range(i);
    p.Fv = @(t) p.m * p.e * p.omega^2 * sin(p.omega * t);
    p.Fh = @(t) p.m * p.e * p.omega^2 * cos(p.omega * t);
    [~, Y] = ode45(@(t,y) equations_of_motion(t,y,p), [0 30], y0, options);
    n_steady = round(5 / 0.005);
    if length(Y) > n_steady
        amp_x_sym(i) = max(Y(end-n_steady:end,1)) - min(Y(end-n_steady:end,1));
    else
        amp_x_sym(i) = NaN;
    end
end

% Asymmetric FRF
p.k1 = @(t) 800; p.k2 = @(t) 500;
for i = 1:length(omega_range)
    p.omega = omega_range(i);
    p.Fv = @(t) p.m * p.e * p.omega^2 * sin(p.omega * t);
    p.Fh = @(t) p.m * p.e * p.omega^2 * cos(p.omega * t);
    [~, Y] = ode45(@(t,y) equations_of_motion(t,y,p), [0 30], y0, options);
    n_steady = round(5 / 0.005);
    if length(Y) > n_steady
        amp_x_asym(i) = max(Y(end-n_steady:end,1)) - min(Y(end-n_steady:end,1));
    else
        amp_x_asym(i) = NaN;
    end
end

figure('Name', 'Frequency Response Comparison', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(omega_range, amp_x_sym, 'k--', 'LineWidth', 1); hold on;
plot(omega_range, amp_x_asym, 'r-', 'LineWidth', 1);
xlabel('\omega (rad/s)', 'FontSize', 10);
ylabel('Amplitude of $x$ (m)', 'Interpreter', 'latex', 'FontSize', 10);
legend('Symmetric', 'Asymmetric', 'Location', 'best');
grid on; box off; set(gca, 'FontSize', 10);

disp('Comparison completed (including accelerations).');