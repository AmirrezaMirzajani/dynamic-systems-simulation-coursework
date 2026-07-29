function dydt = equations_of_motion_timevar(t, y, p)
% State-space form with time-varying m' and k1
% Computes I_m_prime and I_total at each time step

x = y(1);
x_dot = y(2);
theta = y(3);
theta_dot = y(4);

% Evaluate time-dependent parameters
k1_t = p.k1(t);
k2_t = p.k2(t);
m_prime_t = p.m_prime(t);
Fv_t = p.Fv(t);
Fh_t = p.Fh(t);

% Moment of inertia of hollow mass at time t (using current m_prime_t)
A_s = p.L_prime^2;
A_c = pi * p.e^2;
I_m_prime_t = m_prime_t * ( (1/6)*p.L_prime^4 - 0.5*pi*p.e^4 ) / (A_s - A_c);

% Total moment of inertia about G
I_total_t = p.I_beam + I_m_prime_t + m_prime_t * p.d^2;

% Mass matrix (diagonal)
M11 = p.M + m_prime_t;
M22 = I_total_t;

% Damping
C11 = p.c1;

% Stiffness matrix
K11 = k1_t + k2_t;
K12 = (k1_t - k2_t) * (p.L / 2);
K22 = (k1_t + k2_t) * (p.L^2 / 4);

% Generalized forces
Qx = Fv_t;
Qtheta = - p.d * Fh_t;   % negative sign for CCW positive theta

% Accelerations
x_ddot = (Qx - C11*x_dot - K11*x - K12*theta) / M11;
theta_ddot = (Qtheta - K12*x - K22*theta) / M22;

dydt = [x_dot; x_ddot; theta_dot; theta_ddot];
end