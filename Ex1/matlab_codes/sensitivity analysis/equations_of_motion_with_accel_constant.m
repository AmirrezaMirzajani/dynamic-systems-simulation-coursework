function [dydt, x_ddot, theta_ddot] = equations_of_motion_with_accel_constant(t, y, p)
x = y(1);
x_dot = y(2);
theta = y(3);
theta_dot = y(4);

k1 = p.k1;
k2 = p.k2;
m_prime = p.m_prime;
Fv = p.Fv(t);
Fh = p.Fh(t);
I_total = p.I_total;

M11 = p.M + m_prime;
M22 = I_total;
C11 = p.c1;
K11 = k1 + k2;
K12 = (k1 - k2) * (p.L / 2);
K22 = (k1 + k2) * (p.L^2 / 4);
Qx = Fv;
Qtheta = - p.d * Fh;

x_ddot = (Qx - C11*x_dot - K11*x - K12*theta) / M11;
theta_ddot = (Qtheta - K12*x - K22*theta) / M22;

dydt = [x_dot; x_ddot; theta_dot; theta_ddot];
end