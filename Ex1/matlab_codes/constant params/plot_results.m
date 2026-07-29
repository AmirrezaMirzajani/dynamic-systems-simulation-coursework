% plot_results.m
% Plot all required figures: time responses, phase portrait, damper power, and frequency response

% Define black color
myColor = [0 0 0];

% Set default font to Times New Roman
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultLegendFontName', 'Times New Roman');

% Figure size for paper
figWidth = 560;
figHeight = 420;

% ------------------------------------------------------------
% 1. Vertical displacement x(t)
figure('Name', 'x(t)', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T, x, 'Color', myColor, 'LineWidth', 1.5);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$x(t)$ [m]', 'Interpreter', 'latex', 'FontSize', 10);
grid on; box off; set(gca, 'FontSize', 10);

% 2. Vertical velocity x_dot(t)
figure('Name', 'xdot(t)', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T, x_dot, 'Color', myColor, 'LineWidth', 1.5);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$\dot{x}(t)$ [m/s]', 'Interpreter', 'latex', 'FontSize', 10);
grid on; box off; set(gca, 'FontSize', 10);

% 3. Rotation angle theta(t)
figure('Name', 'theta(t)', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T, theta, 'Color', myColor, 'LineWidth', 1.5);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$\theta(t)$ [rad]', 'Interpreter', 'latex', 'FontSize', 10);
grid on; box off; set(gca, 'FontSize', 10);

% 4. Angular velocity theta_dot(t)
figure('Name', 'thetadot(t)', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T, theta_dot, 'Color', myColor, 'LineWidth', 1.5);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$\dot{\theta}(t)$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 10);
grid on; box off; set(gca, 'FontSize', 10);

% 5. Phase portrait (x vs theta)
figure('Name', 'Phase Portrait', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(x, theta, 'Color', myColor, 'LineWidth', 1.5);
xlabel('$x$ [m]', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('$\theta$ [rad]', 'Interpreter', 'latex', 'FontSize', 10);
grid on; box off; set(gca, 'FontSize', 10);


% 7. Frequency response (Amplitude of x vs omega)
% Check if FRF data exist; if not, compute it (or load from file)
if ~exist('omega_range', 'var') || ~exist('amp_x', 'var')
    % Try to load from a saved file (optional)
    if exist('frf_data.mat', 'file')
        load('frf_data.mat', 'omega_range', 'amp_x');
    else
        % Compute FRF on the fly (may take a few seconds)
        disp('Computing frequency response...');
        omega_vec = 1:0.5:50;  % rad/s
        amp_x = zeros(size(omega_vec));
        y0 = [0;0;0;0];
        options = odeset('RelTol',1e-6,'AbsTol',1e-8);
        for i = 1:length(omega_vec)
            p.omega = omega_vec(i);
            p.Fv = @(t) p.m * p.e * p.omega^2 * sin(p.omega * t);
            p.Fh = @(t) p.m * p.e * p.omega^2 * cos(p.omega * t);
            [~, Y] = ode45(@(t,y) equations_of_motion(t,y,p), [0 30], y0, options);
            % Extract steady-state (last 5 seconds)
            fs = 1/0.005;  % sampling frequency (assuming dt=0.005)
            n_steady = round(5 / 0.005);
            if length(Y) > n_steady
                amp_x(i) = max(Y(end-n_steady:end,1)) - min(Y(end-n_steady:end,1));
            else
                amp_x(i) = NaN;
            end
        end
        omega_range = omega_vec;
        % Save for future use
        save('frf_data.mat', 'omega_range', 'amp_x');
    end
end

% Plot FRF if data exists
if exist('omega_range', 'var') && exist('amp_x', 'var') && ~isempty(amp_x)
    figure('Name', 'Frequency Response', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
    plot(omega_range, amp_x, 'Color', myColor, 'LineWidth', 1.5);
    xlabel('\omega (rad/s)', 'FontSize', 10);
    ylabel('Amplitude of $x$ (m)', 'Interpreter', 'latex', 'FontSize', 10);
    grid on; box off; set(gca, 'FontSize', 10);
else
    warning('Could not generate frequency response data. Skipping FRF plot.');
end

% 7. Linear acceleration x_ddot(t)
figure('Name', 'x_ddot(t)', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T, x_ddot, 'Color', myColor, 'LineWidth', 1.5);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$\ddot{x}(t)$ [m/s$^2$]', 'Interpreter', 'latex', 'FontSize', 10);
grid on; box off; set(gca, 'FontSize', 10);

% 8. Angular acceleration theta_ddot(t)
figure('Name', 'theta_ddot(t)', 'NumberTitle', 'off', 'Position', [100 100 figWidth figHeight]);
plot(T, theta_ddot, 'Color', myColor, 'LineWidth', 1.5);
xlabel('Time (s)', 'FontSize', 10);
ylabel('$\ddot{\theta}(t)$ [rad/s$^2$]', 'Interpreter', 'latex', 'FontSize', 10);
grid on; box off; set(gca, 'FontSize', 10);