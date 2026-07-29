%% Natural Frequency Analysis of a Simply Supported Beam with a Discrete Spring
% This script calculates the natural frequencies of a simply supported
% Euler-Bernoulli beam with a discrete spring attached at the beam midpoint.
% The assumed mode method is used with real beam coordinate x.

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
a = 0.25*L;              % Force location [m]

xSpring = L/2;           % Spring location [m]
maxModes = 4;            % Maximum number of assumed modes

%% Preallocation
freqRR = NaN(maxModes,maxModes);
omegaRR = NaN(maxModes,maxModes);
freqError = NaN(maxModes,maxModes);

%% Classical natural frequencies of the beam without spring
modeNumber = (1:maxModes).';

omegaTheory = (modeNumber*pi/L).^2 .* sqrt(E*I/(rho*A));
freqTheory = omegaTheory/(2*pi);

%% Natural frequency calculation using the assumed mode method
for n = 1:maxModes

    [M,K] = buildMK(n,rho,A,E,I,L,k,xSpring);

    [eigVec,eigVal] = eig(K,M);

    lambda = real(diag(eigVal));
    lambda(lambda < 0 & abs(lambda) < 1e-9) = 0;

    [lambda,sortIndex] = sort(lambda,'ascend');
    eigVec = eigVec(:,sortIndex);

    omega = sqrt(lambda);
    freq = omega/(2*pi);

    omegaRR(n,1:n) = omega.';
    freqRR(n,1:n) = freq.';

end

%% Convergence error relative to the four-mode model
freqRef = freqRR(maxModes,:);

for n = 1:maxModes
    for j = 1:n
        freqError(n,j) = abs((freqRR(n,j) - freqRef(j))/freqRef(j))*100;
    end
end

%% Four-mode final matrices and eigenvectors
[M4,K4] = buildMK(maxModes,rho,A,E,I,L,k,xSpring);

[V4,D4] = eig(K4,M4);

lambda4 = real(diag(D4));
lambda4(lambda4 < 0 & abs(lambda4) < 1e-9) = 0;

[lambda4,sortIndex4] = sort(lambda4,'ascend');
V4 = V4(:,sortIndex4);

omegaFinal = sqrt(lambda4);
freqFinal = omegaFinal/(2*pi);

%% Mass normalization of eigenvectors
for r = 1:maxModes
    V4(:,r) = V4(:,r)/sqrt(V4(:,r)'*M4*V4(:,r));
end

%% Print model data
fprintf('\n==================== Model Data ====================\n');
fprintf('rho      = %.2f kg/m^3\n',rho);
fprintf('E        = %.4e Pa\n',E);
fprintf('L        = %.4f m\n',L);
fprintf('b        = %.4f m\n',b);
fprintf('h        = %.4f m\n',h);
fprintf('A        = %.6e m^2\n',A);
fprintf('I        = %.6e m^4\n',I);
fprintf('k        = %.4f N/m\n',k);
fprintf('c        = %.4f N.s/m\n',c);
fprintf('a        = %.4f m\n',a);
fprintf('xSpring  = %.4f m\n',xSpring);

%% Create convergence table
convergenceTable = table((1:maxModes).',freqRR(:,1),freqRR(:,2),freqRR(:,3),freqRR(:,4), ...
    'VariableNames',{'Number_of_Modes','f1_Hz','f2_Hz','f3_Hz','f4_Hz'});

fprintf('\n==================== Natural Frequency Convergence ====================\n');
disp(convergenceTable);

%% Create convergence error table
errorTable = table((1:maxModes).',freqError(:,1),freqError(:,2),freqError(:,3),freqError(:,4), ...
    'VariableNames',{'Number_of_Modes','Error_f1_percent','Error_f2_percent','Error_f3_percent','Error_f4_percent'});

fprintf('\n==================== Convergence Error Relative to n = 4 ====================\n');
disp(errorTable);

%% Create final comparison table
springEffectPercent = 100*(freqFinal - freqTheory)./freqTheory;

finalTable = table(modeNumber,freqTheory,freqFinal,springEffectPercent, ...
    'VariableNames',{'Mode','Classical_Beam_Hz','Beam_with_Spring_Hz','Spring_Effect_percent'});

fprintf('\n==================== Final Four-Mode Frequency Results ====================\n');
disp(finalTable);

%% Save tables
writetable(convergenceTable,'natural_frequency_results.xlsx','Sheet','Convergence');
writetable(errorTable,'natural_frequency_results.xlsx','Sheet','Error');
writetable(finalTable,'natural_frequency_results.xlsx','Sheet','Final_Comparison');

%% Figure 1: Mode shapes of the final four-mode model
x = linspace(0,L,800).';
Phi = zeros(length(x),maxModes);

for i = 1:maxModes
    Phi(:,i) = sin(i*pi*x/L);
end

modeShape = Phi*V4;

for r = 1:maxModes
    [~,idMax] = max(abs(modeShape(:,r)));
    if modeShape(idMax,r) < 0
        modeShape(:,r) = -modeShape(:,r);
    end
    modeShape(:,r) = modeShape(:,r)/max(abs(modeShape(:,r)));
end

figure('Color','w','Name','Mode Shapes');

hold on;
grid on;
box on;

for r = 1:maxModes
    plot(x,modeShape(:,r),'LineWidth',1.8);
end

xline(xSpring,'--','Spring location','LineWidth',1.4);

xlabel('Beam coordinate x [m]');
ylabel('Normalized mode shape');
title('Normalized mode shapes of the four-mode model');
legend('Mode 1','Mode 2','Mode 3','Mode 4','Location','best');

exportgraphics(gcf,'Fig_01_Mode_Shapes.png','Resolution',300);

%% Local function
function [M,K] = buildMK(n,rho,A,E,I,L,k,xSpring)

    M = zeros(n,n);
    K = zeros(n,n);

    for i = 1:n
        for j = 1:n

            if i == j
                M(i,j) = rho*A*L/2;
                K(i,j) = E*I*i^4*pi^4/(2*L^3);
            end

            phi_i_spring = sin(i*pi*xSpring/L);
            phi_j_spring = sin(j*pi*xSpring/L);

            K(i,j) = K(i,j) + k*phi_i_spring*phi_j_spring;

        end
    end

end