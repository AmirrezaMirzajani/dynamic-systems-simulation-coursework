function p = params_timevar()
% Parameters for time-varying simulation (k1(t) and m'(t))
% No precomputation of I_m_prime or I_total because they depend on m'(t)

% Constant parameters
p.M = 2.0;          % mass of beam (kg)
p.I_beam = 0.5;     % moment of inertia of beam about G (kg.m^2)
p.L_prime = 0.2;    % side length of hollow square (m)
p.e = 0.05;         % radius of removed circle & eccentricity of m (m)
p.d = 0.1;          % vertical distance from G to center of m' and rotor (m)
p.L = 0.8;          % beam length (distance between two springs) (m)
p.c1 = 10;          % damper coefficient (N.s/m)
p.m = 0.1;          % unbalanced rotor mass (kg)
p.omega = 20;       % rotor angular speed (rad/s)

% Time-varying parameters (function handles)
% p.k1 = @(t) 800 * (1 - 0.03*t);   % 3% reduction per second, after 10s becomes 560
p.k1 = @(t) 800 * exp(-0.08*t);  
p.k2 = @(t) 800;                       % constant stiffness (N/m)
p.m_prime = @(t) 0.5 + 0.2 * (1./(1+exp(-10*(t-2)))); % smooth increase 0.5->0.7

% Excitation forces (function handles)
p.Fv = @(t) p.m * p.e * p.omega^2 * sin(p.omega * t);
p.Fh = @(t) p.m * p.e * p.omega^2 * cos(p.omega * t);

end