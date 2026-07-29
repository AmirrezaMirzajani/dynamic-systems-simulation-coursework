function p = parameters()
    % System parameters for the mass-pendulum system
    p.M = 2.0;      % kg   upper mass
    p.m = 0.5;      % kg   pendulum mass
    p.L = 0.8;      % m    pendulum length
    p.k = 50;       % N/m  linear spring stiffness
    p.c = 3;        % N.s/m linear damping coefficient
    p.kt = 5;       % N.m/rad torsional spring stiffness
    p.ct = 0.5;     % N.m.s/rad torsional damping coefficient
    p.g = 9.81;     % m/s^2 gravity acceleration
end