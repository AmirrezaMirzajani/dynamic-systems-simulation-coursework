function xdot = dynamics_nonlinear(t, x, p, F_fun, T_fun)
    % State vector: x = [x; x_dot; theta; theta_dot]
    % Returns: xdot = [x_dot; x_ddot; theta_dot; theta_ddot]
    % Nonlinear equations from derived EOM (Eq. 29)
    
    % Extract states
    x_pos = x(1);
    x_vel = x(2);
    theta = x(3);
    theta_vel = x(4);
    
    % External forces evaluated at time t
    F = F_fun(t);
    T = T_fun(t);
    
    % Unpack parameters
    M = p.M;
    m = p.m;
    L = p.L;
    k = p.k;
    c = p.c;
    kt = p.kt;
    ct = p.ct;
    g = p.g;
    
    % Mass matrix M_mat * [x_ddot; theta_ddot] = rhs
    M11 = M + m;
    M12 = 0.5 * m * L * cos(theta);
    M21 = M12;
    M22 = (1/3) * m * L^2;
    
    % Right-hand side vector (including damping, stiffness, gravity, and nonlinear couplings)
    rhs1 = F - c*x_vel - k*x_pos + 0.5*m*L*sin(theta)*theta_vel^2;
    rhs2 = T - ct*theta_vel - kt*theta - 0.5*m*g*L*sin(theta);
    
    % Solve for accelerations
    M_mat = [M11, M12; M21, M22];
    rhs = [rhs1; rhs2];
    acc = M_mat \ rhs;   % [x_ddot; theta_ddot]
    
    x_ddot = acc(1);
    theta_ddot = acc(2);
    
    % Return derivatives
    xdot = [x_vel; x_ddot; theta_vel; theta_ddot];
end