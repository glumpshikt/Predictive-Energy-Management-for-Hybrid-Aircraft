%% Model Predictive Control for Predictive Energy Management of Hybrid Aircraft
clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_mot = 4;                          % number of propulsion systems
m_bat = 8000;                       % battery mass for all systems (kg)
dens = 0.875;                       % battery density (MJ/kg)
SOC_max = dens*m_bat/n_mot*10^6;    % max battery stored energy per system (J)
U_max = 800;                        % battery maximum open circuit voltage (V)
R_max = 3.5e-2;                     % battery eqvl circuit resistance (ohm)
tag = '';                           % file tag
savecsv = false;                    % saves data to csv if true

d = 60;                     % sample time (s)
T = 3600;                   % duration of flight (s)
t=(0:d:T)';                 % time at EM sampling time (s)
N = length(t);              % total number of sampling points

% Generate coefficients of quadratic for R and U
voc = readmatrix('data/Voc_normalised.csv');
R = readmatrix('data/Rbat_normalised.csv');

U_coeff = U_max*polyfit(voc(:,1), voc(:,2), 2);
R_coeff = R_max*polyfit(R(:,1), R(:,2), 2);
clear voc R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generating noise
var_w = 2e-5;                                  % SOC process noise variance
L = 30;                                        % Number of times to repeat the optimisation
K = 20;                                        % Number of scenarios to generate within the optimisation
w = wgn(N, L, var_w, 'linear');                % SOC process noise

var_eps1 = 5000;                               % measurement noise variance for I
var_eps2 = 1000;                               % measurement noise variance for V
eps1 = wgn(N, L, var_eps1, 'linear');          % measurement noise for I
eps2 = wgn(N, L, var_eps2, 'linear');          % measurement noise for V

% Initialise covariance matrices
Ps = zeros(N,L);
Ps(1,:) = var_w*ones(1,L);
Qs = var_w;
Rs = [var_eps1 0; 0 var_eps2];

% Initialise estimates
E = zeros(N,L);
E_true = zeros(N,L);
m = zeros(N, L);
P_c = zeros(N-1, L);
phi = zeros(N-1 ,L);

E(1,:) = 1487.5*10^6/SOC_max*ones(1,L);
E_true(1,:) = E(1,:) + w(N,:); 
m(1,:) = 42000*ones(1,L);
x=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Stepping Through Time %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for l=1:L
    for i=1:N-1
        % generate K noise samples of length N-i+1
        E_k = E(i,l)*ones(K,1) + wgn(K, 1, Ps(i,l), 'linear');
        w_k = wgn(K, N-i+1, var_w, 'linear');

        % get optimum P_b
        if x ~= 0
            x(1:N-i+2:end) = [];
        end
        [P_c(i,l), phi(i,l), x] = ...
            test_IPOPT_MPC(i, E_k*SOC_max/10^6, m(i,l), w_k*SOC_max/10^6, E(i,l)*SOC_max/10^6, x);

        m(i+1,l) = m(i,l) - phi(i,l)*d*n_mot;

        E_true(i+1,l) = stateupdate(E_true(i,l), w(i,l), P_c(i,l), d, SOC_max, U_coeff, R_coeff);

        % State of charge estimation
        E_est = stateupdate(E(i,l), 0, P_c(i,l), d, SOC_max, U_coeff, R_coeff);

        A = calcA(E(i,l), P_c(i,l), d, SOC_max, U_coeff, R_coeff);
        B = calcB(E(i,l), P_c(i,l), d, SOC_max, U_coeff, R_coeff);
        C = calcC(E_est, P_c(i,l), U_coeff, R_coeff);
        D = calcD(E_est, P_c(i,l), U_coeff, R_coeff);

        P_est = A*Ps(i,l)*A+Qs;

        % measurement update
        L = P_est*C.'*inv(C*P_est*C.'+Rs);
        [I_true, V_true] = statemeasure(E_true(i,l), eps1(i,l), eps2(i,l), P_c(i,l), U_coeff, R_coeff);
        y = [I_true; V_true];
        [I_est, V_est] = statemeasure(E_est, 0, 0, P_c(i,l), U_coeff, R_coeff);
        y_est = [I_est; V_est];
        E(i+1,l) = E_est + L*(y - y_est);
        Ps(i+1,l) = P_est - L*C*P_est;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Writing Files to csv %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if savecsv == true 
    filenamem = strcat('data\m_MPC_',tag,'.csv');
    filenameE = strcat('data\E_MPC_',tag,'.csv');
    filenameP_c = strcat('data\P_c_MPC_',tag,'.csv');
    filenamephi = strcat('data\phi_MPC_',tag,'.csv');
    filenameE_true = strcat('data\E_true_MPC_',tag,'.csv');
    
    writematrix(m, filenamem);
    writematrix(E, filenameE);
    writematrix(P_c, filenameP_c);
    writematrix(phi, filenamephi);
    writematrix(E_true, filenameE_true);
end