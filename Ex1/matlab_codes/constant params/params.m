function p = params(varargin)
% Define system parameters (constant or time-varying)
% Use varargin to allow future parameter overrides

% ----- Constant parameters (reasonable assumed values) -----
p.M = 2.0;          % mass of beam (kg)
p.I_beam = 0.5;     % moment of inertia of beam about G (kg.m^2)
p.L_prime = 0.2;    % side length of hollow square (m)
p.e = 0.05;         % radius of removed circle & eccentricity of m (m)
p.d = 0.1;          % vertical distance from G to center of m' and rotor (m)
p.L = 0.8;          % beam length (distance between two springs) (m)
p.c1 = 10;          % damper coefficient (N.s/m)
p.m = 0.1;          % unbalanced rotor mass (kg)
p.omega = 20;       % rotor angular speed (rad/s)

% Spring stiffnesses (can be time-dependent)
const=0;
if const==0
p.k1 = @(t) 800;    % left spring (N/m)
p.k2 = @(t) 800;    % right spring (N/m)
p.m_prime = 0.5;    % hollow mass (kg)
elseif const==1   
p.k1 = @(t) 800 + 80*sin(2*pi*1*t);   % oscillating stiffness (N/m)
p.k2 = @(t) 800;                       % constant stiffness (N/m)
p.m_prime = @(t) 0.5 + 0.2 * (1./(1+exp(-10*(t-2)))); % smooth increase 0.5->0.7 around t=2s
end


% ----- Moment of inertia of hollow mass (square with circular hole) -----
A_s = p.L_prime^2;
A_c = pi * p.e^2;
if A_s <= A_c
    error('Circle area larger than square. Reduce e.');
end
p.I_m_prime = p.m_prime * ( (1/6)*p.L_prime^4 - 0.5*pi*p.e^4 ) / (A_s - A_c);

% ----- Total moment of inertia about G (beam + m' shifted by distance d) -----
p.I_total = p.I_beam + p.I_m_prime + p.m_prime * p.d^2;

% ----- Excitation forces (vertical and horizontal) -----
p.Fv = @(t) p.m * p.e * p.omega^2 * sin(p.omega * t);
p.Fh = @(t) p.m * p.e * p.omega^2 * cos(p.omega * t);

% Overwrite fields if provided via varargin (optional)
for i = 1:2:nargin
    if ischar(varargin{i})
        p.(varargin{i}) = varargin{i+1};
    end
end
end