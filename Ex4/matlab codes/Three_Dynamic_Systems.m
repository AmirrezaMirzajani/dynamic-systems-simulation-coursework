%% Random Gaussian White Noise Analysis for Three Dynamic Systems
% This script applies Gaussian white noise excitation to three dynamic systems,
% solves the time-domain response, and performs frequency analysis using
% FFT and Welch PSD estimation.

clear;
clc;
close all;

%% General simulation settings
rng(25);

dt = 0.001;                  % Time step [s]
tEnd = 20.0;                 % Simulation time [s]
t = (0:dt:tEnd).';           % Time vector [s]
fs = 1/dt;                   % Sampling frequency [Hz]
N = length(t);

fMinPeak = 0.2;              % Minimum frequency for peak search [Hz]
fMaxPeak = 250;              % Maximum frequency for peak search [Hz]
fMaxPlot = 250;              % Maximum frequency for plots [Hz]
nPeaks = 6;                  % Number of dominant frequency peaks
minPeakDistanceHz = 1.0;     % Minimum distance between selected peaks [Hz]

outputExcel = 'random_noise_three_systems_results.xlsx';

if isfile(outputExcel)
    delete(outputExcel);
end

%% ========================================================================
%  System 1: Rigid beam-spring-damper system with translational and rotational DOF
%  Generalized coordinates: x(t), theta(t)
%  Random inputs: vertical white-noise force and horizontal white-noise force
% ========================================================================

%% System 1 parameters
p1.M = 2.0;                  % Beam mass [kg]
p1.I = 0.5;                  % Beam mass moment of inertia about G [kg.m^2]
p1.mp = 0.5;                 % Hollow attached mass [kg]
p1.Lp = 0.2;                 % Side length of hollow square mass [m]
p1.e = 0.05;                 % Hole radius and rotor eccentricity [m]
p1.d = 0.1;                  % Vertical distance from G [m]
p1.L = 0.8;                  % Beam length between springs [m]
p1.c1 = 10;                  % Central viscous damping coefficient [N.s/m]

p1.k1 = 800;                 % Left spring stiffness [N/m]
p1.k2 = 800;                 % Right spring stiffness [N/m]

p1.Im = p1.mp*((1/6)*p1.Lp^4 - 0.5*pi*p1.e^4)/(p1.Lp^2 - pi*p1.e^2);
p1.Itotal = p1.I + p1.Im + p1.mp*p1.d^2;

sigmaFv1 = 10;               % Standard deviation of vertical white-noise force [N]
sigmaFh1 = 10;               % Standard deviation of horizontal white-noise force [N]

Fv1 = sigmaFv1*randn(N,1);
Fh1 = sigmaFh1*randn(N,1);

z01 = zeros(4,1);

[~,z1] = ode45(@(time,z) system1ODE(time,z,p1,t,Fv1,Fh1),t,z01);

x1 = z1(:,1);
xd1 = z1(:,2);
theta1 = z1(:,3);
thetad1 = z1(:,4);

freqNat1 = naturalFreqSystem1(p1);

fprintf('\n==================== System 1 Natural Frequencies ====================\n');
disp(freqNat1);

%% System 1 frequency analysis
x1_mm = x1*1000;
theta1_deg = theta1*180/pi;

[fftF_x1,fftA_x1] = singleSidedFFT(x1_mm,fs);
[fftF_th1,fftA_th1] = singleSidedFFT(theta1_deg,fs);

[psdF_x1,psdA_x1] = welchPSD(x1_mm,fs);
[psdF_th1,psdA_th1] = welchPSD(theta1_deg,fs);

peaksFFT_x1 = createPeakTable("System 1","x","FFT amplitude [mm]",fftF_x1,fftA_x1,nPeaks,fMinPeak,fMaxPeak,minPeakDistanceHz);
peaksFFT_th1 = createPeakTable("System 1","theta","FFT amplitude [deg]",fftF_th1,fftA_th1,nPeaks,fMinPeak,fMaxPeak,minPeakDistanceHz);

peaksPSD_x1 = createPeakTable("System 1","x","PSD [mm^2/Hz]",psdF_x1,psdA_x1,nPeaks,fMinPeak,fMaxPeak,minPeakDistanceHz);
peaksPSD_th1 = createPeakTable("System 1","theta","PSD [deg^2/Hz]",psdF_th1,psdA_th1,nPeaks,fMinPeak,fMaxPeak,minPeakDistanceHz);

summary1 = createSummaryTable("System 1", ...
    ["x";"theta"], ...
    ["mm";"deg"], ...
    {x1_mm,theta1_deg}, ...
    {peaksFFT_x1,peaksFFT_th1}, ...
    {peaksPSD_x1,peaksPSD_th1});

fprintf('\n==================== System 1 Response Summary ====================\n');
disp(summary1);

fprintf('\n==================== System 1 FFT Peaks ====================\n');
disp([peaksFFT_x1;peaksFFT_th1]);

fprintf('\n==================== System 1 PSD Peaks ====================\n');
disp([peaksPSD_x1;peaksPSD_th1]);

%% System 1 plots
figure('Color','w','Name','System 1 Inputs');

tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(t,Fv1,'LineWidth',1.0);
grid on;
box on;
xlabel('Time [s]');
ylabel('F_v [N]');
title('System 1 vertical Gaussian white-noise force');

nexttile;
plot(t,Fh1,'LineWidth',1.0);
grid on;
box on;
xlabel('Time [s]');
ylabel('F_h [N]');
title('System 1 horizontal Gaussian white-noise force');

exportgraphics(gcf,'Fig_01_System1_Noise_Inputs.png','Resolution',300);

figure('Color','w','Name','System 1 Time Response');

tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(t,x1_mm,'LineWidth',1.2);
grid on;
box on;
xlabel('Time [s]');
ylabel('x [mm]');
title('System 1 displacement response');

nexttile;
plot(t,theta1_deg,'LineWidth',1.2);
grid on;
box on;
xlabel('Time [s]');
ylabel('\theta [deg]');
title('System 1 rotational response');

exportgraphics(gcf,'Fig_02_System1_Time_Response.png','Resolution',300);

figure('Color','w','Name','System 1 Frequency Analysis');

tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(fftF_x1,fftA_x1,'LineWidth',1.2);
grid on;
box on;
xlim([0 fMaxPlot]);
xlabel('Frequency [Hz]');
ylabel('Amplitude [mm]');
title('System 1 FFT of x');

nexttile;
plot(fftF_th1,fftA_th1,'LineWidth',1.2);
grid on;
box on;
xlim([0 fMaxPlot]);
xlabel('Frequency [Hz]');
ylabel('Amplitude [deg]');
title('System 1 FFT of \theta');

nexttile;
plot(psdF_x1,psdA_x1,'LineWidth',1.2);
grid on;
box on;
xlim([0 fMaxPlot]);
xlabel('Frequency [Hz]');
ylabel('PSD [mm^2/Hz]');
title('System 1 Welch PSD of x');

nexttile;
plot(psdF_th1,psdA_th1,'LineWidth',1.2);
grid on;
box on;
xlim([0 fMaxPlot]);
xlabel('Frequency [Hz]');
ylabel('PSD [deg^2/Hz]');
title('System 1 Welch PSD of \theta');

exportgraphics(gcf,'Fig_03_System1_Frequency_Analysis.png','Resolution',300);

%% ========================================================================
%  System 2: Nonlinear cart-pendulum-like mass-rod system
%  Generalized coordinates: x(t), theta(t)
%  Random inputs: horizontal white-noise force and optional torque noise
% ========================================================================

%% System 2 parameters
p2.M = 2.0;                  % Cart mass [kg]
p2.m = 0.5;                  % Rod mass [kg]
p2.L = 0.8;                  % Rod length [m]
p2.k = 50;                   % Linear spring stiffness [N/m]
p2.c = 3;                    % Linear damping coefficient [N.s/m]
p2.kt = 5;                   % Torsional spring stiffness [N.m/rad]
p2.ct = 0.5;                 % Torsional damping coefficient [N.m.s/rad]
p2.g = 9.81;                 % Gravity acceleration [m/s^2]
p2.I = (1/3)*p2.m*p2.L^2;    % Rod inertia about hinge [kg.m^2]

sigmaF2 = 2.0;               % Standard deviation of horizontal white-noise force [N]
sigmaT2 = 0.5;               % Standard deviation of white-noise torque [N.m]
applyTorqueNoise2 = true;    % Set false if only force input is required

F2 = sigmaF2*randn(N,1);

if applyTorqueNoise2
    T2 = sigmaT2*randn(N,1);
else
    T2 = zeros(N,1);
end

z02 = zeros(4,1);

[~,z2] = ode45(@(time,z) system2ODE(time,z,p2,t,F2,T2),t,z02);

x2 = z2(:,1);
xd2 = z2(:,2);
theta2 = z2(:,3);
thetad2 = z2(:,4);

freqNat2 = naturalFreqSystem2(p2);

fprintf('\n==================== System 2 Linearized Natural Frequencies ====================\n');
disp(freqNat2);

%% System 2 frequency analysis
x2_mm = x2*1000;
theta2_deg = theta2*180/pi;

[fftF_x2,fftA_x2] = singleSidedFFT(x2_mm,fs);
[fftF_th2,fftA_th2] = singleSidedFFT(theta2_deg,fs);

[psdF_x2,psdA_x2] = welchPSD(x2_mm,fs);
[psdF_th2,psdA_th2] = welchPSD(theta2_deg,fs);

peaksFFT_x2 = createPeakTable("System 2","x","FFT amplitude [mm]",fftF_x2,fftA_x2,nPeaks,fMinPeak,fMaxPeak,minPeakDistanceHz);
peaksFFT_th2 = createPeakTable("System 2","theta","FFT amplitude [deg]",fftF_th2,fftA_th2,nPeaks,fMinPeak,fMaxPeak,minPeakDistanceHz);

peaksPSD_x2 = createPeakTable("System 2","x","PSD [mm^2/Hz]",psdF_x2,psdA_x2,nPeaks,fMinPeak,fMaxPeak,minPeakDistanceHz);
peaksPSD_th2 = createPeakTable("System 2","theta","PSD [deg^2/Hz]",psdF_th2,psdA_th2,nPeaks,fMinPeak,fMaxPeak,minPeakDistanceHz);

summary2 = createSummaryTable("System 2", ...
    ["x";"theta"], ...
    ["mm";"deg"], ...
    {x2_mm,theta2_deg}, ...
    {peaksFFT_x2,peaksFFT_th2}, ...
    {peaksPSD_x2,peaksPSD_th2});

fprintf('\n==================== System 2 Response Summary ====================\n');
disp(summary2);

fprintf('\n==================== System 2 FFT Peaks ====================\n');
disp([peaksFFT_x2;peaksFFT_th2]);

fprintf('\n==================== System 2 PSD Peaks ====================\n');
disp([peaksPSD_x2;peaksPSD_th2]);

%% System 2 plots
figure('Color','w','Name','System 2 Inputs');

tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(t,F2,'LineWidth',1.0);
grid on;
box on;
xlabel('Time [s]');
ylabel('F(t) [N]');
title('System 2 Gaussian white-noise force');

nexttile;
plot(t,T2,'LineWidth',1.0);
grid on;
box on;
xlabel('Time [s]');
ylabel('T(t) [N.m]');
title('System 2 Gaussian white-noise torque');

exportgraphics(gcf,'Fig_04_System2_Noise_Inputs.png','Resolution',300);

figure('Color','w','Name','System 2 Time Response');

tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(t,x2_mm,'LineWidth',1.2);
grid on;
box on;
xlabel('Time [s]');
ylabel('x [mm]');
title('System 2 displacement response');

nexttile;
plot(t,theta2_deg,'LineWidth',1.2);
grid on;
box on;
xlabel('Time [s]');
ylabel('\theta [deg]');
title('System 2 angular response');

exportgraphics(gcf,'Fig_05_System2_Time_Response.png','Resolution',300);

figure('Color','w','Name','System 2 Frequency Analysis');

tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(fftF_x2,fftA_x2,'LineWidth',1.2);
grid on;
box on;
xlim([0 fMaxPlot]);
xlabel('Frequency [Hz]');
ylabel('Amplitude [mm]');
title('System 2 FFT of x');

nexttile;
plot(fftF_th2,fftA_th2,'LineWidth',1.2);
grid on;
box on;
xlim([0 fMaxPlot]);
xlabel('Frequency [Hz]');
ylabel('Amplitude [deg]');
title('System 2 FFT of \theta');

nexttile;
plot(psdF_x2,psdA_x2,'LineWidth',1.2);
grid on;
box on;
xlim([0 fMaxPlot]);
xlabel('Frequency [Hz]');
ylabel('PSD [mm^2/Hz]');
title('System 2 Welch PSD of x');

nexttile;
plot(psdF_th2,psdA_th2,'LineWidth',1.2);
grid on;
box on;
xlim([0 fMaxPlot]);
xlabel('Frequency [Hz]');
ylabel('PSD [deg^2/Hz]');
title('System 2 Welch PSD of \theta');

exportgraphics(gcf,'Fig_06_System2_Frequency_Analysis.png','Resolution',300);

%% ========================================================================
%  System 3: Simply supported Euler-Bernoulli beam with discrete spring and damper
%  Modal model with four assumed modes
%  Random input: point Gaussian white-noise force
% ========================================================================

%% System 3 parameters
p3.rho = 7850;               % Density [kg/m^3]
p3.E = 210e9;                % Young's modulus [Pa]
p3.L = 1.0;                  % Beam length [m]
p3.b = 0.030;                % Beam width [m]
p3.h = 0.005;                % Beam height [m]
p3.A = p3.b*p3.h;            % Cross-sectional area [m^2]
p3.I = p3.b*p3.h^3/12;       % Second moment of area [m^4]
p3.k = 5000;                 % Discrete spring stiffness [N/m]
p3.c = 5;                    % Discrete damping coefficient [N.s/m]
p3.a = 0.25*p3.L;            % Force location from left support [m]
p3.xSpring = p3.L/2;         % Spring location [m]
p3.xDamper = p3.L/2;         % Damper location [m]
p3.xMid = p3.L/2;            % Midpoint observation location [m]
p3.xForce = p3.a;            % Force-point observation location [m]
p3.nModes = 4;               % Number of assumed modes

sigmaF3 = 1.0;               % Standard deviation of point white-noise force [N]
F3 = sigmaF3*randn(N,1);

[M3,C3,K3] = buildBeamMCK(p3);
B3 = modalForceVector(p3);

z03 = zeros(2*p3.nModes,1);

[~,z3] = ode45(@(time,z) beamModalODE(time,z,M3,C3,K3,B3,t,F3),t,z03);

q3 = z3(:,1:p3.nModes);

w3Mid = reconstructBeamDisplacement(q3,p3.nModes,p3.xMid,p3.L);
w3Force = reconstructBeamDisplacement(q3,p3.nModes,p3.xForce,p3.L);

freqNat3 = naturalFrequenciesFromMK(M3,K3);

fprintf('\n==================== System 3 Natural Frequencies ====================\n');
disp(freqNat3);

fprintf('\n==================== System 3 Modal Force Participation Vector ====================\n');
disp(B3);

%% System 3 frequency analysis
w3Mid_mm = w3Mid*1000;
w3Force_mm = w3Force*1000;

[fftF_mid3,fftA_mid3] = singleSidedFFT(w3Mid_mm,fs);
[fftF_force3,fftA_force3] = singleSidedFFT(w3Force_mm,fs);

[psdF_mid3,psdA_mid3] = welchPSD(w3Mid_mm,fs);
[psdF_force3,psdA_force3] = welchPSD(w3Force_mm,fs);

peaksFFT_mid3 = createPeakTable("System 3","w at L/2","FFT amplitude [mm]",fftF_mid3,fftA_mid3,nPeaks,fMinPeak,fMaxPeak,minPeakDistanceHz);
peaksFFT_force3 = createPeakTable("System 3","w at a","FFT amplitude [mm]",fftF_force3,fftA_force3,nPeaks,fMinPeak,fMaxPeak,minPeakDistanceHz);

peaksPSD_mid3 = createPeakTable("System 3","w at L/2","PSD [mm^2/Hz]",psdF_mid3,psdA_mid3,nPeaks,fMinPeak,fMaxPeak,minPeakDistanceHz);
peaksPSD_force3 = createPeakTable("System 3","w at a","PSD [mm^2/Hz]",psdF_force3,psdA_force3,nPeaks,fMinPeak,fMaxPeak,minPeakDistanceHz);

summary3 = createSummaryTable("System 3", ...
    ["w at L/2";"w at a"], ...
    ["mm";"mm"], ...
    {w3Mid_mm,w3Force_mm}, ...
    {peaksFFT_mid3,peaksFFT_force3}, ...
    {peaksPSD_mid3,peaksPSD_force3});

fprintf('\n==================== System 3 Response Summary ====================\n');
disp(summary3);

fprintf('\n==================== System 3 FFT Peaks ====================\n');
disp([peaksFFT_mid3;peaksFFT_force3]);

fprintf('\n==================== System 3 PSD Peaks ====================\n');
disp([peaksPSD_mid3;peaksPSD_force3]);

%% System 3 plots
figure('Color','w','Name','System 3 Input');

plot(t,F3,'LineWidth',1.0);
grid on;
box on;
xlabel('Time [s]');
ylabel('F(t) [N]');
title('System 3 Gaussian white-noise point force');

exportgraphics(gcf,'Fig_07_System3_Noise_Input.png','Resolution',300);

figure('Color','w','Name','System 3 Time Response');

tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(t,w3Mid_mm,'LineWidth',1.2);
grid on;
box on;
xlabel('Time [s]');
ylabel('w(L/2,t) [mm]');
title('System 3 midpoint displacement response');

nexttile;
plot(t,w3Force_mm,'LineWidth',1.2);
grid on;
box on;
xlabel('Time [s]');
ylabel('w(a,t) [mm]');
title('System 3 force-point displacement response');

exportgraphics(gcf,'Fig_08_System3_Time_Response.png','Resolution',300);

figure('Color','w','Name','System 3 Frequency Analysis');

tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(fftF_mid3,fftA_mid3,'LineWidth',1.2);
grid on;
box on;
xlim([0 fMaxPlot]);
xlabel('Frequency [Hz]');
ylabel('Amplitude [mm]');
title('System 3 FFT of w at L/2');

nexttile;
plot(fftF_force3,fftA_force3,'LineWidth',1.2);
grid on;
box on;
xlim([0 fMaxPlot]);
xlabel('Frequency [Hz]');
ylabel('Amplitude [mm]');
title('System 3 FFT of w at a');

nexttile;
plot(psdF_mid3,psdA_mid3,'LineWidth',1.2);
grid on;
box on;
xlim([0 fMaxPlot]);
xlabel('Frequency [Hz]');
ylabel('PSD [mm^2/Hz]');
title('System 3 Welch PSD of w at L/2');

nexttile;
plot(psdF_force3,psdA_force3,'LineWidth',1.2);
grid on;
box on;
xlim([0 fMaxPlot]);
xlabel('Frequency [Hz]');
ylabel('PSD [mm^2/Hz]');
title('System 3 Welch PSD of w at a');

exportgraphics(gcf,'Fig_09_System3_Frequency_Analysis.png','Resolution',300);

%% ========================================================================
%  Save all result tables
% ========================================================================

allSummary = [summary1;summary2;summary3];

allFFTPeaks = [peaksFFT_x1;peaksFFT_th1; ...
               peaksFFT_x2;peaksFFT_th2; ...
               peaksFFT_mid3;peaksFFT_force3];

allPSDPeaks = [peaksPSD_x1;peaksPSD_th1; ...
               peaksPSD_x2;peaksPSD_th2; ...
               peaksPSD_mid3;peaksPSD_force3];

writetable(allSummary,outputExcel,'Sheet','Summary');
writetable(allFFTPeaks,outputExcel,'Sheet','FFT_Peaks');
writetable(allPSDPeaks,outputExcel,'Sheet','PSD_Peaks');

writetable(table(freqNat1(:),'VariableNames',{'Frequency_Hz'}),outputExcel,'Sheet','System1_NaturalFreq');
writetable(table(freqNat2(:),'VariableNames',{'Frequency_Hz'}),outputExcel,'Sheet','System2_NaturalFreq');
writetable(table(freqNat3(:),'VariableNames',{'Frequency_Hz'}),outputExcel,'Sheet','System3_NaturalFreq');

fprintf('\n==================== All Summary Results ====================\n');
disp(allSummary);

fprintf('\nResults saved to: %s\n',outputExcel);

%% ========================================================================
%  Local functions
% ========================================================================

function dzdt = system1ODE(time,z,p,tNoise,FvNoise,FhNoise)

    x = z(1);
    xd = z(2);
    theta = z(3);
    thetad = z(4);

    Fv = interp1(tNoise,FvNoise,time,'linear',0);
    Fh = interp1(tNoise,FhNoise,time,'linear',0);

    xdd = (Fv - p.c1*xd - (p.k1+p.k2)*x + (p.L/2)*(p.k1-p.k2)*theta)/(p.M+p.mp);

    thetadd = (p.d*Fh - (p.L/2)*(p.k1-p.k2)*x - (p.L^2/4)*(p.k1+p.k2)*theta)/p.Itotal;

    dzdt = [xd; xdd; thetad; thetadd];

end

function freqHz = naturalFreqSystem1(p)

    Mmat = [p.M+p.mp, 0; ...
            0, p.Itotal];

    Kmat = [p.k1+p.k2, -(p.L/2)*(p.k1-p.k2); ...
            -(p.L/2)*(p.k1-p.k2), (p.L^2/4)*(p.k1+p.k2)];

    eigValue = eig(Kmat,Mmat);
    eigValue = real(eigValue);
    eigValue(eigValue < 0 & abs(eigValue) < 1e-9) = 0;

    freqHz = sqrt(sort(eigValue,'ascend'))/(2*pi);

end

function dzdt = system2ODE(time,z,p,tNoise,FNoise,TNoise)

    x = z(1);
    xd = z(2);
    theta = z(3);
    thetad = z(4);

    F = interp1(tNoise,FNoise,time,'linear',0);
    T = interp1(tNoise,TNoise,time,'linear',0);

    M11 = p.M + p.m;
    M12 = 0.5*p.m*p.L*cos(theta);
    M21 = M12;
    M22 = (1/3)*p.m*p.L^2;

    Mmat = [M11, M12; M21, M22];

    rhs1 = F - p.c*xd - p.k*x + 0.5*p.m*p.L*sin(theta)*thetad^2;
    rhs2 = T - p.ct*thetad - p.kt*theta - 0.5*p.m*p.g*p.L*sin(theta);

    qdd = Mmat \ [rhs1; rhs2];

    xdd = qdd(1);
    thetadd = qdd(2);

    dzdt = [xd; xdd; thetad; thetadd];

end

function freqHz = naturalFreqSystem2(p)

    Mmat = [p.M+p.m, 0.5*p.m*p.L; ...
            0.5*p.m*p.L, (1/3)*p.m*p.L^2];

    Kmat = [p.k, 0; ...
            0, p.kt + 0.5*p.m*p.g*p.L];

    eigValue = eig(Kmat,Mmat);
    eigValue = real(eigValue);
    eigValue(eigValue < 0 & abs(eigValue) < 1e-9) = 0;

    freqHz = sqrt(sort(eigValue,'ascend'))/(2*pi);

end

function [M,C,K] = buildBeamMCK(p)

    n = p.nModes;

    M = zeros(n,n);
    C = zeros(n,n);
    K = zeros(n,n);

    for i = 1:n
        for j = 1:n

            if i == j
                M(i,j) = p.rho*p.A*p.L/2;
                K(i,j) = p.E*p.I*i^4*pi^4/(2*p.L^3);
            end

            phi_i_spring = sin(i*pi*p.xSpring/p.L);
            phi_j_spring = sin(j*pi*p.xSpring/p.L);

            phi_i_damper = sin(i*pi*p.xDamper/p.L);
            phi_j_damper = sin(j*pi*p.xDamper/p.L);

            K(i,j) = K(i,j) + p.k*phi_i_spring*phi_j_spring;
            C(i,j) = C(i,j) + p.c*phi_i_damper*phi_j_damper;

        end
    end

end

function B = modalForceVector(p)

    n = p.nModes;
    B = zeros(n,1);

    for i = 1:n
        B(i) = sin(i*pi*p.a/p.L);
    end

end

function dzdt = beamModalODE(time,z,M,C,K,B,tForce,Fforce)

    n = size(M,1);

    q = z(1:n);
    qd = z(n+1:2*n);

    Fcurrent = interp1(tForce,Fforce,time,'linear',0);

    qdd = M \ (B*Fcurrent - C*qd - K*q);

    dzdt = [qd; qdd];

end

function w = reconstructBeamDisplacement(q,n,xObs,L)

    phiObs = zeros(n,1);

    for i = 1:n
        phiObs(i) = sin(i*pi*xObs/L);
    end

    w = q*phiObs;

end

function freqHz = naturalFrequenciesFromMK(M,K)

    eigValue = eig(K,M);
    eigValue = real(eigValue);
    eigValue(eigValue < 0 & abs(eigValue) < 1e-9) = 0;

    freqHz = sqrt(sort(eigValue,'ascend'))/(2*pi);

end

function [fAxis,ampSpectrum] = singleSidedFFT(signal,fs)

    signal = signal(:);
    signal = signal - mean(signal);

    N = length(signal);

    window = hannWindow(N);
    coherentGain = mean(window);

    Y = fft(signal.*window);

    P2 = abs(Y/(N*coherentGain));
    P1 = P2(1:floor(N/2)+1);

    if length(P1) > 2
        P1(2:end-1) = 2*P1(2:end-1);
    end

    fAxis = fs*(0:floor(N/2))/N;
    ampSpectrum = P1;

end

function [fAxis,PSD] = welchPSD(signal,fs)

    signal = signal(:);
    signal = signal - mean(signal);

    N = length(signal);

    segmentLength = min(round(4*fs),N);
    segmentLength = max(segmentLength,256);

    overlapLength = round(0.5*segmentLength);
    stepLength = segmentLength - overlapLength;

    nfft = 2^nextpow2(segmentLength);

    window = hannWindow(segmentLength);
    U = sum(window.^2);

    startIndex = 1;
    counter = 0;
    PSDacc = zeros(floor(nfft/2)+1,1);

    while (startIndex + segmentLength - 1) <= N

        segment = signal(startIndex:startIndex+segmentLength-1);
        segment = segment - mean(segment);

        X = fft(segment.*window,nfft);

        P2 = (abs(X).^2)/(fs*U);
        P1 = P2(1:floor(nfft/2)+1);

        if length(P1) > 2
            P1(2:end-1) = 2*P1(2:end-1);
        end

        PSDacc = PSDacc + P1;

        counter = counter + 1;
        startIndex = startIndex + stepLength;

    end

    if counter == 0
        PSD = zeros(floor(nfft/2)+1,1);
    else
        PSD = PSDacc/counter;
    end

    fAxis = fs*(0:floor(nfft/2))/nfft;

end

function window = hannWindow(N)

    if N <= 1
        window = ones(N,1);
    else
        n = (0:N-1).';
        window = 0.5*(1 - cos(2*pi*n/(N-1)));
    end

end

function peakTable = createPeakTable(systemName,signalName,amplitudeName,fAxis,responseSpectrum,nPeaks,fMin,fMax,minDistanceHz)

    [peakFreq,peakValue] = topSpectrumPeaks(fAxis,responseSpectrum,nPeaks,fMin,fMax,minDistanceHz);

    n = length(peakFreq);

    peakTable = table( ...
        repmat(string(systemName),n,1), ...
        repmat(string(signalName),n,1), ...
        (1:n).', ...
        peakFreq(:), ...
        peakValue(:), ...
        repmat(string(amplitudeName),n,1), ...
        'VariableNames',{'System','Signal','Peak_Number','Frequency_Hz','Value','Value_Type'});

end

function [peakFreq,peakValue] = topSpectrumPeaks(fAxis,responseSpectrum,nPeaks,fMin,fMax,minDistanceHz)

    fAxis = fAxis(:);
    responseSpectrum = responseSpectrum(:);

    validID = (fAxis >= fMin) & (fAxis <= fMax);

    f = fAxis(validID);
    y = responseSpectrum(validID);

    if isempty(f) || isempty(y)
        peakFreq = NaN;
        peakValue = NaN;
        return;
    end

    if length(y) < 3
        [peakValue,idx] = max(y);
        peakFreq = f(idx);
        return;
    end

    localPeakID = false(size(y));
    localPeakID(2:end-1) = (y(2:end-1) > y(1:end-2)) & (y(2:end-1) >= y(3:end));

    candidateFreq = f(localPeakID);
    candidateValue = y(localPeakID);

    if isempty(candidateValue)
        [candidateValue,idx] = max(y);
        candidateFreq = f(idx);
    end

    [candidateValue,sortID] = sort(candidateValue,'descend');
    candidateFreq = candidateFreq(sortID);

    selectedFreq = [];
    selectedValue = [];

    for i = 1:length(candidateValue)

        if isempty(selectedFreq)
            selectedFreq = candidateFreq(i);
            selectedValue = candidateValue(i);
        else
            distanceOK = all(abs(candidateFreq(i) - selectedFreq) >= minDistanceHz);

            if distanceOK
                selectedFreq = [selectedFreq; candidateFreq(i)];
                selectedValue = [selectedValue; candidateValue(i)];
            end
        end

        if length(selectedFreq) >= nPeaks
            break;
        end

    end

    peakFreq = selectedFreq;
    peakValue = selectedValue;

end

function summaryTable = createSummaryTable(systemName,signalNames,units,signals,fftPeaks,psdPeaks)

    n = length(signals);

    systemColumn = repmat(string(systemName),n,1);
    signalColumn = strings(n,1);
    unitColumn = strings(n,1);

    maxAbsColumn = zeros(n,1);
    rmsColumn = zeros(n,1);
    stdColumn = zeros(n,1);
    meanColumn = zeros(n,1);
    finalAbsColumn = zeros(n,1);
    dominantFFTColumn = zeros(n,1);
    dominantPSDColumn = zeros(n,1);

    for i = 1:n

        signal = signals{i};
        signalColumn(i) = signalNames(i);
        unitColumn(i) = units(i);

        maxAbsColumn(i) = max(abs(signal));
        rmsColumn(i) = sqrt(mean(signal.^2));
        stdColumn(i) = std(signal);
        meanColumn(i) = mean(signal);
        finalAbsColumn(i) = abs(signal(end));

        dominantFFTColumn(i) = fftPeaks{i}.Frequency_Hz(1);
        dominantPSDColumn(i) = psdPeaks{i}.Frequency_Hz(1);

    end

    summaryTable = table(systemColumn,signalColumn,unitColumn,maxAbsColumn,rmsColumn,stdColumn, ...
        meanColumn,finalAbsColumn,dominantFFTColumn,dominantPSDColumn, ...
        'VariableNames',{'System','Signal','Unit','MaxAbs','RMS','STD','Mean','FinalAbs', ...
        'Dominant_FFT_Frequency_Hz','Dominant_PSD_Frequency_Hz'});

end