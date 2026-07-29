%% Free Vibration Analysis of a Simply Supported Beam with Spring and Damper
% This script solves the free vibration response of a simply supported
% Euler-Bernoulli beam with a discrete spring and viscous damper.
% The response is computed for different assumed mode numbers, stiffness
% values, damping values, and FFT analysis.

clear;
clc;
close all;

%% Physical and geometrical properties
rho = 7850;              % Density [kg/m^3]
E   = 210e9;             % Young's modulus [Pa]
L   = 1.0;               % Beam length [m]

b = 0.030;               % Beam width [m]
h = 0.005;               % Beam height [m]

A = b*h;                 % Cross-sectional area [m^2]
I = b*h^3/12;            % Second moment of area [m^4]

k0 = 5000;               % Reference spring stiffness [N/m]
c0 = 5;                  % Reference damping coefficient [N.s/m]

a = 0.25*L;              % Force location [m]
xSpring = L/2;           % Spring location [m]
xDamper = L/2;           % Damper location [m]

xObs = L/2;              % Observation point [m]
maxModes = 4;            % Maximum number of assumed modes

%% Time settings
dt = 0.001;              % Time step [s]
tEnd = 1.5;               % Simulation time [s]
t = 0:dt:tEnd;           % Time vector [s]
fs = 1/dt;               % Sampling frequency [Hz]

%% Initial conditions
qInitialAmplitude = 1e-3;    % Initial modal displacement amplitude [m]

%% Case 1: Free vibration response for n = 1 to n = 4
responseModes = zeros(length(t),maxModes);
summaryModes = zeros(maxModes,3);

for n = 1:maxModes

    [M,C,K] = buildMCK(n,rho,A,E,I,L,k0,c0,xSpring,xDamper);

    q0 = zeros(n,1);
    qd0 = zeros(n,1);

    q0(1) = qInitialAmplitude;

    z0 = [q0; qd0];

    [~,z] = ode45(@(time,z) freeSystemODE(time,z,M,C,K),t,z0);

    q = z(:,1:n);

    wObs = reconstructDisplacement(q,n,xObs,L);

    responseModes(:,n) = wObs;

    summaryModes(n,1) = n;
    summaryModes(n,2) = max(abs(wObs));
    summaryModes(n,3) = rms(wObs);

end

modeSummaryTable = array2table(summaryModes, ...
    'VariableNames',{'Number_of_Modes','MaxAbsDisplacement_m','RMSDisplacement_m'});

fprintf('\n==================== Free Vibration Summary for Different Mode Numbers ====================\n');
disp(modeSummaryTable);

%% Plot 1: Free vibration response for different mode numbers
figure('Color','w','Name','Free Vibration - Mode Number Effect');

hold on;
grid on;
box on;

for n = 1:maxModes
    plot(t,responseModes(:,n),'LineWidth',1.5);
end

xlabel('Time [s]');
ylabel('Displacement at x = L/2 [m]');
title('Free vibration response for different assumed mode numbers');
legend('n = 1','n = 2','n = 3','n = 4','Location','best');

exportgraphics(gcf,'Fig_02_Free_Response_Mode_Number.png','Resolution',300);

%% Case 2: Effect of spring stiffness
kFactors = [0.5 1.0 2.0];
kValues = kFactors*k0;

responseK = zeros(length(t),length(kValues));
summaryK = zeros(length(kValues),5);

for r = 1:length(kValues)

    kCurrent = kValues(r);

    [M,C,K] = buildMCK(maxModes,rho,A,E,I,L,kCurrent,c0,xSpring,xDamper);

    q0 = zeros(maxModes,1);
    qd0 = zeros(maxModes,1);

    q0(1) = qInitialAmplitude;

    z0 = [q0; qd0];

    [~,z] = ode45(@(time,z) freeSystemODE(time,z,M,C,K),t,z0);

    q = z(:,1:maxModes);

    wObs = reconstructDisplacement(q,maxModes,xObs,L);

    responseK(:,r) = wObs;

    freqCurrent = naturalFrequencies(M,K);

    summaryK(r,1) = kCurrent;
    summaryK(r,2) = freqCurrent(1);
    summaryK(r,3) = max(abs(wObs));
    summaryK(r,4) = rms(wObs);
    summaryK(r,5) = abs(wObs(end));

end

kSummaryTable = array2table(summaryK, ...
    'VariableNames',{'k_N_per_m','f1_Hz','MaxAbsDisplacement_m','RMSDisplacement_m','FinalAbsDisplacement_m'});

fprintf('\n==================== Effect of Spring Stiffness ====================\n');
disp(kSummaryTable);

%% Plot 2: Effect of spring stiffness
figure('Color','w','Name','Free Vibration - Spring Stiffness Effect');

hold on;
grid on;
box on;

for r = 1:length(kValues)
    plot(t,responseK(:,r),'LineWidth',1.5);
end

xlabel('Time [s]');
ylabel('Displacement at x = L/2 [m]');
title('Effect of spring stiffness on free vibration response');
legend('k = 0.5 k_0','k = k_0','k = 2 k_0','Location','best');

exportgraphics(gcf,'Fig_03_Free_Response_k_Effect.png','Resolution',300);

%% Case 3: Effect of damping coefficient
cFactors = [0.5 1.0 2.0];
cValues = cFactors*c0;

responseC = zeros(length(t),length(cValues));
summaryC = zeros(length(cValues),4);

for r = 1:length(cValues)

    cCurrent = cValues(r);

    [M,C,K] = buildMCK(maxModes,rho,A,E,I,L,k0,cCurrent,xSpring,xDamper);

    q0 = zeros(maxModes,1);
    qd0 = zeros(maxModes,1);

    q0(1) = qInitialAmplitude;

    z0 = [q0; qd0];

    [~,z] = ode45(@(time,z) freeSystemODE(time,z,M,C,K),t,z0);

    q = z(:,1:maxModes);

    wObs = reconstructDisplacement(q,maxModes,xObs,L);

    responseC(:,r) = wObs;

    summaryC(r,1) = cCurrent;
    summaryC(r,2) = max(abs(wObs));
    summaryC(r,3) = rms(wObs);
    summaryC(r,4) = abs(wObs(end));

end

cSummaryTable = array2table(summaryC, ...
    'VariableNames',{'c_Ns_per_m','MaxAbsDisplacement_m','RMSDisplacement_m','FinalAbsDisplacement_m'});

fprintf('\n==================== Effect of Damping Coefficient ====================\n');
disp(cSummaryTable);

%% Plot 3: Effect of damping coefficient
figure('Color','w','Name','Free Vibration - Damping Effect');

hold on;
grid on;
box on;

for r = 1:length(cValues)
    plot(t,responseC(:,r),'LineWidth',1.5);
end

xlabel('Time [s]');
ylabel('Displacement at x = L/2 [m]');
title('Effect of damping coefficient on free vibration response');
legend('c = 0.5 c_0','c = c_0','c = 2 c_0','Location','best');

exportgraphics(gcf,'Fig_04_Free_Response_c_Effect.png','Resolution',300);

%% Case 4: FFT of the reference free vibration response
wReference = responseModes(:,maxModes);

[fAxis,ampSpectrum] = singleSidedFFT(wReference,fs);

fftTable = table(fAxis(:),ampSpectrum(:), ...
    'VariableNames',{'Frequency_Hz','Amplitude_m'});

fprintf('\n==================== Dominant FFT Peaks of Free Vibration Response ====================\n');

[peakValues,peakLocations] = findpeaks(ampSpectrum,fAxis, ...
    'SortStr','descend', ...
    'NPeaks',5);

peakTable = table(peakLocations(:),peakValues(:), ...
    'VariableNames',{'Frequency_Hz','Amplitude_m'});

disp(peakTable);

%% Plot 4: FFT of free vibration response
figure('Color','w','Name','FFT - Free Vibration Response');

plot(fAxis,ampSpectrum,'LineWidth',1.5);
grid on;
box on;

xlabel('Frequency [Hz]');
ylabel('Amplitude [m]');
title('FFT of free vibration response at x = L/2');
xlim([0 250]);

exportgraphics(gcf,'Fig_05_FFT_Free_Response.png','Resolution',300);

%% Save tables
writetable(modeSummaryTable,'free_vibration_results.xlsx','Sheet','Mode_Number');
writetable(kSummaryTable,'free_vibration_results.xlsx','Sheet','Spring_Effect');
writetable(cSummaryTable,'free_vibration_results.xlsx','Sheet','Damping_Effect');
writetable(peakTable,'free_vibration_results.xlsx','Sheet','FFT_Peaks');

%% Local functions
function [M,C,K] = buildMCK(n,rho,A,E,I,L,k,c,xSpring,xDamper)

    M = zeros(n,n);
    C = zeros(n,n);
    K = zeros(n,n);

    for i = 1:n
        for j = 1:n

            if i == j
                M(i,j) = rho*A*L/2;
                K(i,j) = E*I*i^4*pi^4/(2*L^3);
            end

            phi_i_spring = sin(i*pi*xSpring/L);
            phi_j_spring = sin(j*pi*xSpring/L);

            phi_i_damper = sin(i*pi*xDamper/L);
            phi_j_damper = sin(j*pi*xDamper/L);

            K(i,j) = K(i,j) + k*phi_i_spring*phi_j_spring;
            C(i,j) = C(i,j) + c*phi_i_damper*phi_j_damper;

        end
    end

end

function dzdt = freeSystemODE(~,z,M,C,K)

    n = size(M,1);

    q = z(1:n);
    qd = z(n+1:2*n);

    qdd = M \ (-C*qd - K*q);

    dzdt = [qd; qdd];

end

function w = reconstructDisplacement(q,n,xObs,L)

    phiObs = zeros(n,1);

    for i = 1:n
        phiObs(i) = sin(i*pi*xObs/L);
    end

    w = q*phiObs;

end

function freqHz = naturalFrequencies(M,K)

    eigValue = eig(K,M);
    eigValue = real(eigValue);
    eigValue(eigValue < 0 & abs(eigValue) < 1e-9) = 0;

    omega = sqrt(sort(eigValue,'ascend'));
    freqHz = omega/(2*pi);

end

function [fAxis,ampSpectrum] = singleSidedFFT(signal,fs)

    signal = signal(:);
    signal = signal - mean(signal);

    N = length(signal);

    Y = fft(signal);

    P2 = abs(Y/N);
    P1 = P2(1:floor(N/2)+1);

    P1(2:end-1) = 2*P1(2:end-1);

    fAxis = fs*(0:floor(N/2))/N;
    ampSpectrum = P1;

end