function [dydt, x_ddot, theta_ddot] = equations_of_motion_with_accel_timevar(t, y, p)
% Compute dydt and accelerations for time-varying system

x = y(1);
x_dot = y(2);
theta = y(3);
theta_dot = y(4);

k1_t = p.k1(t);
k2_t = p.k2(t);
m_prime_t = p.m_prime(t);
Fv_t = p.Fv(t);
Fh_t = p.Fh(t);

A_s = p.L_prime^2;
A_c = pi * p.e^2;
I_m_prime_t = m_prime_t * ( (1/6)*p.L_prime^4 - 0.5*pi*p.e^4 ) / (A_s - A_c);
I_total_t = p.I_beam + I_m_prime_t + m_prime_t * p.d^2;

M11 = p.M + m_prime_t;
M22 = I_total_t;
C11 = p.c1;
K11 = k1_t + k2_t;
K12 = (k1_t - k2_t) * (p.L / 2);
K22 = (k1_t + k2_t) * (p.L^2 / 4);
Qx = Fv_t;
Qtheta = - p.d * Fh_t;

x_ddot = (Qx - C11*x_dot - K11*x - K12*theta) / M11;
theta_ddot = (Qtheta - K12*x - K22*theta) / M22;

dydt = [x_dot; x_ddot; theta_dot; theta_ddot];
end