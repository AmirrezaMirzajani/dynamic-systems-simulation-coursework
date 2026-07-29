clear; close all; clc;

%% Load parameters
p = parameters();

%% Force functions (both zero for case 1)
F_fun = @(t) 0;
T_fun = @(t) 0;

%% Three initial condition scenarios: [x0; xd0; theta0; thetad0]
IC1 = [0.1; 0; 0.2; 0];   % Scenario 1
IC2 = [0.2; 0; 0.5; 0];   % Scenario 2
IC3 = [-0.1; 0; 0.7; 0];  % Scenario 3

ICs = {IC1, IC2, IC3};
scenarioNames = {'Scenario 1', 'Scenario 2', 'Scenario 3'};

%% Time span
tspan = [0, 10];

%% Solve ODE for each scenario
sols = cell(1,3);
for i = 1:3
    [t, y] = ode45(@(t,y) dynamics_nonlinear(t,y,p,F_fun,T_fun), tspan, ICs{i});
    sols{i}.t = t;
    sols{i}.x = y(:,1);
    sols{i}.theta = y(:,3);
end

%% Figure 1: Displacement of mass M (x vs t) - all scenarios
figure;
hold on;
for i = 1:3
    plot(sols{i}.t, sols{i}.x, 'LineWidth', 1.5);
end
hold off;
grid on;
xlabel('Time (s)');
ylabel('x(t) (m)');
legend(scenarioNames, 'Location', 'best');
% No title here (you will add manually in Persian)

%% Figure 2: Angle of pendulum (theta vs t) - all scenarios
figure;
hold on;
for i = 1:3
    plot(sols{i}.t, sols{i}.theta, 'LineWidth', 1.5);
end
hold off;
grid on;
xlabel('Time (s)');
ylabel('\theta(t) (rad)');
legend(scenarioNames, 'Location', 'best');
% No title here

%% Figure 3: Phase portrait for Scenario 3 (x vs dx/dt and theta vs dtheta/dt)
% We need velocity data from solver - re-run for scenario 3 or store earlier
% I'll recompute scenario 3 with full state (or we already have in sols but we only stored x,theta)
% Better to store full state in sols from beginning. Let me correct:

% Clear and re-solve properly storing full state
clear sols;
for i = 1:3
    [t, y] = ode45(@(t,y) dynamics_nonlinear(t,y,p,F_fun,T_fun), tspan, ICs{i});
    sols{i}.t = t;
    sols{i}.x = y(:,1);
    sols{i}.xd = y(:,2);
    sols{i}.theta = y(:,3);
    sols{i}.thetad = y(:,4);
end

% Now phase portrait for scenario 3
figure;
plot(sols{3}.x, sols{3}.xd, 'b-', 'LineWidth', 1.2); hold on;
plot(sols{3}.theta, sols{3}.thetad, 'r-', 'LineWidth', 1.2);
hold off;
grid on;
xlabel('Displacement / Angle');
ylabel('Velocity');
legend('x vs dx/dt', '\theta vs d\theta/dt', 'Location', 'best');
% No title here