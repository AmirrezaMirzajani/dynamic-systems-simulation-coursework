%% Forced Vibration Analysis of a Simply Supported Beam with Spring and Damper
% This script solves the forced vibration response of a simply supported
% Euler-Bernoulli beam with a discrete spring and viscous damper.
% Different excitation cases are considered, including pulse, harmonic,
% multi-sine, and Gaussian white noise forces.

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

k = 5000;                % Spring stiffness [N/m]
c = 5;                   % Damping coefficient [N.s/m]

a = 0.25*L;              % Force location from the left support [m]
xSpring = L/2;           % Spring location [m]
xDamper = L/2;           % Damper location [m]

xMid = L/2;              % Midpoint observation location [m]
xForce = a;              % Force-point observation location [m]

nModes = 4;              % Number of assumed modes

%% Time settings
dt = 0.001;              % Time step [s]
tEnd = 1.0;              % Simulation time [s]
t = (0:dt:tEnd).';       % Time vector [s]
fs = 1/dt;               % Sampling frequency [Hz]

%% Build system matrices
[M,C,K] = buildMCK(nModes,rho,A,E,I,L,k,c,xSpring,xDamper);

freqNat = naturalFrequencies(M,K);
omegaNat = 2*pi*freqNat;

fprintf('\n==================== Natural Frequencies Used in Forced Vibration ====================\n');
for i = 1:nModes
    fprintf('Mode %d: %.4f Hz\n',i,freqNat(i));
end

%% Generalized force vector
B = zeros(nModes,1);

for i = 1:nModes
    B(i) = sin(i*pi*a/L);
end

fprintf('\n==================== Modal Force Participation Vector ====================\n');
disp(B);

%% Initial conditions
q0 = zeros(nModes,1);
qd0 = zeros(nModes,1);
z0 = [q0; qd0];

%% Force case 1: Pulse force
F0Pulse = 10;            % Pulse amplitude [N]
tp = 0.02;               % Pulse duration [s]

F_pulse = zeros(size(t));
F_pulse(t <= tp) = F0Pulse;

%% Force case 2: Harmonic force at the first natural frequency
F0Sin = 1.0;             % Harmonic force amplitude [N]
fSin = freqNat(1);       % Excitation frequency [Hz]
omegaSin = 2*pi*fSin;    % Excitation angular frequency [rad/s]

F_sine = F0Sin*sin(omegaSin*t);

%% Force case 3: Multi-sine pseudo-random force
rng(10);

multiAmp = [1.0 0.8 0.6 0.4];                      % Force amplitudes [N]
multiFreq = [7.0 freqNat(1) 55.0 freqNat(3)];       % Force frequencies [Hz]
multiPhase = 2*pi*rand(1,length(multiAmp));         % Phase angles [rad]

F_multisine = zeros(size(t));

for r = 1:length(multiAmp)
    F_multisine = F_multisine + multiAmp(r)*sin(2*pi*multiFreq(r)*t + multiPhase(r));
end

%% Force case 4: Gaussian white noise force
rng(20);

sigmaNoise = 1.0;       % Standard deviation of white noise force [N]
F_noise = sigmaNoise*randn(size(t));

%% Store main force cases
caseNames = {'Pulse';'Sine at f1';'Multi-sine';'White noise'};
forceCases = {F_pulse,F_sine,F_multisine,F_noise};

nCases = length(forceCases);

wMidCases = zeros(length(t),nCases);
wForceCases = zeros(length(t),nCases);

summaryMain = zeros(nCases,7);

%% Solve main forced vibration cases
for r = 1:nCases

    Fcurrent = forceCases{r};

    [~,z] = ode45(@(time,z) forcedSystemODE(time,z,M,C,K,B,t,Fcurrent),t,z0);

    q = z(:,1:nModes);

    wMid = reconstructDisplacement(q,nModes,xMid,L);
    wForcePoint = reconstructDisplacement(q,nModes,xForce,L);

    wMidCases(:,r) = wMid;
    wForceCases(:,r) = wForcePoint;

    [fAxis,ampMid] = singleSidedFFT(wMid,fs);
    [peakFreq,peakAmp] = topSpectrumPeaks(fAxis,ampMid,1,0.5,300);

    summaryMain(r,1) = r;
    summaryMain(r,2) = max(abs(Fcurrent));
    summaryMain(r,3) = max(abs(wMid));
    summaryMain(r,4) = sqrt(mean(wMid.^2));
    summaryMain(r,5) = max(abs(wForcePoint));
    summaryMain(r,6) = sqrt(mean(wForcePoint.^2));
    summaryMain(r,7) = peakFreq(1);

end

mainSummaryTable = array2table(summaryMain, ...
    'VariableNames',{'Case_Number','Max_Force_N','Max_Mid_Displacement_m','RMS_Mid_Displacement_m', ...
    'Max_ForcePoint_Displacement_m','RMS_ForcePoint_Displacement_m','Dominant_Response_Frequency_Hz'});

mainSummaryTable.Case_Name = caseNames;
mainSummaryTable = movevars(mainSummaryTable,'Case_Name','After','Case_Number');

fprintf('\n==================== Forced Vibration Summary for Main Excitation Cases ====================\n');
disp(mainSummaryTable);

%% Force case 5: Harmonic excitation frequency variation
sweepFreq = [0.5*freqNat(1), freqNat(1), 1.5*freqNat(1), freqNat(2)];
sweepLabel = {'0.5 f1';'f1';'1.5 f1';'f2'};

nSweep = length(sweepFreq);

wMidSweep = zeros(length(t),nSweep);
wForceSweep = zeros(length(t),nSweep);

summarySweep = zeros(nSweep,8);

for r = 1:nSweep

    fExc = sweepFreq(r);
    omegaExc = 2*pi*fExc;

    F_sweep = F0Sin*sin(omegaExc*t);

    [~,z] = ode45(@(time,z) forcedSystemODE(time,z,M,C,K,B,t,F_sweep),t,z0);

    q = z(:,1:nModes);

    wMid = reconstructDisplacement(q,nModes,xMid,L);
    wForcePoint = reconstructDisplacement(q,nModes,xForce,L);

    wMidSweep(:,r) = wMid;
    wForceSweep(:,r) = wForcePoint;

    steadyStart = round(0.60*length(t));

    summarySweep(r,1) = r;
    summarySweep(r,2) = fExc;
    summarySweep(r,3) = fExc/freqNat(1);
    summarySweep(r,4) = max(abs(wMid));
    summarySweep(r,5) = max(abs(wMid(steadyStart:end)));
    summarySweep(r,6) = sqrt(mean(wMid.^2));
    summarySweep(r,7) = max(abs(wForcePoint));
    summarySweep(r,8) = max(abs(wForcePoint(steadyStart:end)));

end

sweepSummaryTable = array2table(summarySweep, ...
    'VariableNames',{'Case_Number','Excitation_Frequency_Hz','Frequency_Ratio_to_f1', ...
    'Max_Mid_Displacement_m','Steady_Max_Mid_Displacement_m','RMS_Mid_Displacement_m', ...
    'Max_ForcePoint_Displacement_m','Steady_Max_ForcePoint_Displacement_m'});

sweepSummaryTable.Excitation_Case = sweepLabel;
sweepSummaryTable = movevars(sweepSummaryTable,'Excitation_Case','After','Case_Number');

fprintf('\n==================== Harmonic Excitation Frequency Sweep ====================\n');
disp(sweepSummaryTable);

%% FFT peak table for main forced cases
peakData = [];

for r = 1:nCases

    [fAxis,ampMid] = singleSidedFFT(wMidCases(:,r),fs);
    [peakFreq,peakAmp] = topSpectrumPeaks(fAxis,ampMid,3,0.5,300);

    for p = 1:length(peakFreq)
        peakData = [peakData; r, p, peakFreq(p), peakAmp(p)];
    end

end

fftPeakTable = array2table(peakData, ...
    'VariableNames',{'Case_Number','Peak_Number','Frequency_Hz','Amplitude_m'});

caseNameColumn = strings(height(fftPeakTable),1);

for r = 1:height(fftPeakTable)
    caseNameColumn(r) = caseNames{fftPeakTable.Case_Number(r)};
end

fftPeakTable.Case_Name = caseNameColumn;
fftPeakTable = movevars(fftPeakTable,'Case_Name','After','Case_Number');

fprintf('\n==================== FFT Peaks for Main Forced Vibration Cases ====================\n');
disp(fftPeakTable);

%% Figure 1: Input force histories
figure('Color','w','Name','Input Force Histories');

tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

for r = 1:nCases
    nexttile;
    plot(t,forceCases{r},'LineWidth',1.2);
    grid on;
    box on;
    ylabel('F(t) [N]');
    title(caseNames{r});
end

xlabel('Time [s]');

exportgraphics(gcf,'Fig_06_Input_Force_Histories.png','Resolution',300);

%% Figure 2: Midpoint responses for main forced cases
figure('Color','w','Name','Forced Responses at Midpoint');

hold on;
grid on;
box on;

for r = 1:nCases
    plot(t,wMidCases(:,r),'LineWidth',1.3);
end

xlabel('Time [s]');
ylabel('Displacement at x = L/2 [m]');
title('Forced vibration response at beam midpoint');
legend(caseNames,'Location','best');

exportgraphics(gcf,'Fig_07_Forced_Response_Midpoint.png','Resolution',300);

%% Figure 3: FFT of midpoint responses
figure('Color','w','Name','FFT of Forced Responses');

hold on;
grid on;
box on;

for r = 1:nCases
    [fAxis,ampMid] = singleSidedFFT(wMidCases(:,r),fs);
    plot(fAxis,ampMid,'LineWidth',1.3);
end

xlabel('Frequency [Hz]');
ylabel('Amplitude [m]');
title('FFT of forced vibration responses at beam midpoint');
legend(caseNames,'Location','best');
xlim([0 250]);

exportgraphics(gcf,'Fig_08_FFT_Forced_Response_Midpoint.png','Resolution',300);

%% Figure 4: Harmonic sweep time responses at midpoint
figure('Color','w','Name','Harmonic Sweep Time Responses');

tiledlayout(nSweep,1,'TileSpacing','compact','Padding','compact');

for r = 1:nSweep

    nexttile;

    plot(t,wMidSweep(:,r)*1000,'LineWidth',1.2);
    grid on;
    box on;

    ylabel('w [mm]');
    title(['Excitation frequency: ', sweepLabel{r}]);

end

xlabel('Time [s]');

exportgraphics(gcf,'Fig_09_Harmonic_Sweep_Time_Responses.png','Resolution',300);

%% Save result tables
writetable(mainSummaryTable,'forced_vibration_results.xlsx','Sheet','Main_Cases');
writetable(sweepSummaryTable,'forced_vibration_results.xlsx','Sheet','Frequency_Sweep');
writetable(fftPeakTable,'forced_vibration_results.xlsx','Sheet','FFT_Peaks');

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

function dzdt = forcedSystemODE(time,z,M,C,K,B,tForce,Fforce)

    n = size(M,1);

    q = z(1:n);
    qd = z(n+1:2*n);

    Fcurrent = interp1(tForce,Fforce,time,'linear','extrap');

    qdd = M \ (B*Fcurrent - C*qd - K*q);

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

function [peakFreq,peakAmp] = topSpectrumPeaks(fAxis,ampSpectrum,nPeaks,fMin,fMax)

    fAxis = fAxis(:);
    ampSpectrum = ampSpectrum(:);

    validID = (fAxis >= fMin) & (fAxis <= fMax);

    f = fAxis(validID);
    a = ampSpectrum(validID);

    if length(a) < 3
        peakFreq = f;
        peakAmp = a;
        return;
    end

    localPeakID = false(size(a));
    localPeakID(2:end-1) = (a(2:end-1) > a(1:end-2)) & (a(2:end-1) >= a(3:end));

    peakFreqAll = f(localPeakID);
    peakAmpAll = a(localPeakID);

    if isempty(peakAmpAll)
        [peakAmpAll,idx] = max(a);
        peakFreqAll = f(idx);
    end

    [peakAmpSorted,sortID] = sort(peakAmpAll,'descend');
    peakFreqSorted = peakFreqAll(sortID);

    nReturn = min(nPeaks,length(peakAmpSorted));

    peakFreq = peakFreqSorted(1:nReturn);
    peakAmp = peakAmpSorted(1:nReturn);

end