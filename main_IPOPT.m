%% Energy management for hybrid electric aircraft
% Optimise the power split between the gas turbine and
% electric motor of a parallel hybrid propulsion system to meet the drive
% power demand for a given flight profile. Open loop implementation of the
% energy management algorithm in the paper "Predictive energy management
% for hybrid electric aircraft propulsion systems".

clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solver parameters
open_circuit_deg = 2;       % 0: Constant, >0: number of deg of polynomial
resistance_deg = 2;         % 0: Constant, >0: number of deg of polynomial
tag = '';                   % file tag
savecsv = false;            % saves data to csv if true
d = 24;                     % sample time (s)
T = 3600;                   % duration of flight (s)

% Physical parameters
n_mot = 4;                  % number of propulsion systems
m_bat = 8000;               % battery mass for all systems (kg)
dens = 0.875;               % battery density (MJ/kg)
SOC_max = dens*m_bat/n_mot; % max battery stored energy per system (MJ)
R_max = 3.5e-2;             % battery eqvl circuit resistance (ohm)
U_max = 800;                % battery maximum open circuit voltage (V)
E_low = 0.2*SOC_max;        % lower bound on battery energy per system (MJ)
E_up = 0.85*SOC_max;        % upper bound on battery energy per system (MJ)
P_em_low = 0;               % lower bound on e.m. power per shaft (MW)
P_em_up = 5;                % upper bound on e.m. power per shaft (MW)
P_gt_low = 0;               % lower bound on gas turbine power (MW)
P_gt_up = 5;                % upper bound on gas turbine power (MW)
compCl = [0.11 0.43];       % lift coefficients: [b0(-); b1(deg^-1)]
compCd = [0.0282 0.0042...
    5.3613e-04];            % drag coeff.: [a0(-); a1(deg^-1); a2(deg^-2)]
initTAS = 131.9263;         % initial true airspeed (m/s)
g = 9.81;                   % acceleration due to gravity (m/s^2)
S = 77.3;                   % characteristic surface area (m^2)
M = 42000;                  % aircraft maximum take-off weight (kg)
eta_w = 0.5;                % battery charge up efficiency (-)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Trajectory planning %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Flight parameters
vtas_up = 190;                           % cruise TAS (m/s)
h_up = 12750;                            % cruise altitude (m)

% Flight profile generation
vp=[initTAS  vtas_up...
    vtas_up initTAS]';                   % TAS data points (m/s)
tvp = [0 500 2000 T]';                   % TAS time data points (s)
hp = [0 h_up h_up 0];                    % height data points (m)
thp = [0 1000 2000 T]';                  % height time data points (s)
t=(0:d:T)';                              % time at EM sampling time (s)
v = interp1(tvp, vp, t, 'pchip');        % interp. of TAS profile (m/s)
height = interp1(thp, hp, t, 'pchip');   % interp. of height profile (m)
h_dot=diff(height)./diff(t);             % vertical speed (m/s)
gamma=asin(h_dot./v(1:end-1));
gamma=[gamma; gamma(end)];               % flight path angle (rad) (level)
v_dot = diff(v)./diff(t);
v_dot = [v_dot; v_dot(end)];             % acceleration profile (m/s2)
gamma_dot = diff(gamma)./diff(t);
gamma_dot = [gamma_dot; gamma_dot(end)]; % rate of flight path angle (rad/s)

rho_air = get_rho(height);               % air density at ref alt. (kg/m^3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Power profile & losses %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Coefficients for drive power
N = length(t);
Phi = ones(length(t), 1);

% Get drive power coeff. in P_drv = eta(:,1).*m.^2 + eta(:,2).*m + eta(:,3)
[eta2, eta1, eta0] = f_drv(v, v_dot, gamma, gamma_dot, compCd, compCl,...
    g, S, rho_air);                     % (MW per shaft)
eta = [eta2, eta1, eta0];

% gt fuel coeff. in ph = beta1*Peng + beta0 (1st order fit)
beta = [0.0821 0.0327].*ones(N, 2);     % (beta = [beta1; beta0])

% e.m. loss map coeff. in Pc = alpha2*P_gt^2 + alpha1*P_gt + alpha0
alpha = [0 1 0].*ones(N, 3);            % (alpha = [alpha2; alpha1; alpha0])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Power constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Bounds on variables
% Bounds on electric motor input (electrical) power P_c
P_c_low = h(P_em_low, alpha);
P_c_up = h(P_em_up, alpha);

% Bounds on fuel consumption
phi_low = f(P_gt_low*Phi, beta);
phi_up = f(P_gt_up*Phi, beta);

% Estimate bounds for P_b
P_b_low_est = inv_f_P_b(P_c_low, U_max, R_max);
P_b_up_est1 = inv_f_P_b(P_c_up, U_max, R_max);
P_b_up_est2 = U_max^2/(2*R_max)*10^-6;
P_b_up_est = min(P_b_up_est1, P_b_up_est2*Phi);

%% scale factors for parameters
E_scale = E_up;
U_scale = U_max;
P_b_scale = max(P_b_up_est);
phi_scale = max(phi_up);
M_scale = M;
R_scale = R_max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Battery Coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

voc = readmatrix('data/Voc_normalised.csv');
R = readmatrix('data/Rbat_normalised.csv');

% get polynomial coefficients for v_oc
if open_circuit_deg == 0
    v_coeff = 1;
else
    v_coeff = polyfit(voc(:,1), voc(:,2), open_circuit_deg);
end

% get polynomial coefficients for R
if resistance_deg == 0
    r_coeff = 1;
else
    r_coeff = polyfit(R(:,1), R(:,2), open_circuit_deg);
end

clear voc R;

%% Initialise values for IPOPT
E_0 = E_up*Phi/E_scale;
U_0 = U_max*Phi/U_scale;
m_0 = M*Phi/M_scale;
R_0 = R_max*Phi/R_scale;

x0 = [E_0; U_0; P_b_up_est/P_b_scale; phi_up/phi_scale; m_0; R_0];

auxdata = {N, alpha, beta, eta, U_max, R_max, SOC_max, ...         %(1-7)
    M, E_up, P_em_up, v_coeff, r_coeff, d, n_mot, ...              %(8-14)
    E_scale, U_scale, P_b_scale, phi_scale, M_scale, R_scale, ...  %(15-20)
    P_em_low};                                                     %(21)

%% we take the optimisation variable as x = [ E, U, P_b, phi, m, R ];
% The callback functions.
funcs.objective         = @objective;
funcs.constraints       = @constraints;
funcs.gradient          = @gradipopt;
funcs.jacobian          = @jacobian;
funcs.jacobianstructure = @jacobiansparse;

% Lower bound on the variables.
options.lb = [E_low*ones(N,1)/E_scale; -Inf*ones(N,1); ...
    -Inf*ones(N,1); phi_low/phi_scale; -Inf*ones(N, 1); -Inf*ones(N,1)];
% Upper bound on the variables.
options.ub = [E_up*ones(N,1)/E_scale; Inf*ones(N,1); ...
    Inf*ones(N,1); phi_up/phi_scale; Inf*ones(N,1); Inf*ones(N,1)];
% Lower bounds on the constraint functions.
options.cl = [zeros(N,1); zeros(N,1); zeros(N,1); ...
    zeros(N,1); -Inf*ones(N,1); -Inf*ones(N,1); zeros(N,1); zeros(N,1)];
% Upper bounds on the constraint functions.
options.cu = [Inf*ones(N,1); zeros(N,1); zeros(N,1); ...
    zeros(N,1); zeros(N,1); zeros(N,1); zeros(N,1); Inf*ones(N,1)];

% Set the IPOPT options.
options.ipopt.mu_strategy = 'adaptive';
options.ipopt.tol         = 1e-10;
options.ipopt.print_level = 5;
options.ipopt.hessian_approximation = 'limited-memory';
options.auxdata = {N, alpha, beta, eta, U_max, R_max, SOC_max, ... %(1-7)
    M, E_up, P_em_up, v_coeff, r_coeff, d, n_mot, ...              %(8-14)
    E_scale, U_scale, P_b_scale, phi_scale, M_scale, R_scale, ...  %(15-20)
    P_em_low};

%% Run IPOPT.
[x, info] = ipopt_auxdata(x0,funcs,options);

%% Write data for IPOPT
E = x(1:N)*E_scale;
U = x(N+1: 2*N)*U_scale;
P_b = [x(2*N+1: 3*N-1); x(3*N-1)]*P_b_scale;
phi = [x(3*N+1: 4*N); x(4*N-1)]*phi_scale;
m = x(4*N+1: 5*N)*M_scale;
R = x(5*N+1: 6*N)*R_scale;

if savecsv == true 
    filenameE = strcat('data\E_',tag,'.csv');
    filenameU = strcat('data\U_',tag,'.csv');
    filenameP_b = strcat('data\P_b_',tag,'.csv');
    filenamephi = strcat('data\phi_',tag,'.csv');
    filenamem = strcat('data\m_',tag,'.csv');
    filenameR = strcat('data\R_',tag,'.csv');

    writematrix(E, filenameE);
    writematrix(U, filenameU);
    writematrix(P_b, filenameP_b);
    writematrix(phi, filenamephi);
    writematrix(m, filenamem);
    writematrix(R, filenameR);
end

%% Calculate power split
P_em = inv_h(inv_f_P_b(P_b, U, R), alpha);

% Hybridisation ratio
P_drv = f_P_drv(m, eta);
H_ratio = P_em./P_drv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
burnt_fuel = M - m(end);
fprintf('Amount of fuel burnt = %f kg', burnt_fuel);

%% plotting figures
figure
hold on
plot(t, P_drv)
plot(t, P_em)
legend('$P_{drv}$', '$P_{em}$','interpreter', 'latex')
xlim([0, 3600])
ylabel('Power (MW)', 'interpreter', 'latex')
xlabel('time (s)', 'interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex')
grid on
hold off

function o = objective(x, auxdata)
N = cell2mat(auxdata(1));
M = cell2mat(auxdata(8));
M_scale = cell2mat(auxdata(19));
o = M/M_scale-x(N*5);
end

function g = gradipopt(~, auxdata)
N = cell2mat(auxdata(1));
g = [zeros(N*5-1,1); -1; zeros(N,1)];
end
