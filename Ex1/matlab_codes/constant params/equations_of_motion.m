function dydt = equations_of_motion(t, y, p)
% State-space form of equations of motion
% y = [x; x_dot; theta; theta_dot]
% dydt = [x_dot; x_ddot; theta_dot; theta_ddot]

x = y(1);
x_dot = y(2);
theta = y(3);
theta_dot = y(4);

% Time-varying parameters
k1_t = p.k1(t);
k2_t = p.k2(t);
Fv_t = p.Fv(t);
Fh_t = p.Fh(t);

% Mass matrix (diagonal)
M11 = p.M + p.m_prime;
M22 = p.I_total;

% Damping (only in x)
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