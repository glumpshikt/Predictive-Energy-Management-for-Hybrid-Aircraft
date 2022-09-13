function [P_c, phi, x] = test_IPOPT_MPC(timestep, E_init, M_init, w_k, E_est, x0)
%% Solver parameters
open_circuit_deg = 2;       % 0: Constant, >0: number of deg of polynomial
resistance_deg = 2;         % 0: Constant, >0: number of deg of polynomial
d = 60;                     % sample time (s)
T = 3600;                   % duration of flight (s)
K = length(E_init);         % number of estimates

% Physical parameters
n_mot = 4;                  % number of propulsion systems
SOC_max = 1750;             % Maximum state of charge (MJ)
E_low = 0.2*SOC_max;        % lower bound on battery energy per system (MJ)
E_up = max(0.85*SOC_max, max(E_init)-min(min(w_k)));        % upper bound on battery energy per system (MJ)
R_max = 0.035;              % Battery internal resistance (Ohms)
U_max = 800;                % Maximum open-circuit voltage (V)
P_em_low = 0;               % lower bound on e.m. power per shaft (MW)
P_em_up = 5;                % upper bound on e.m. power per shaft (MW)
P_gt_low = 0;               % lower bound on gas turbine power (MW)
P_gt_up = 5;                % upper bound on gas turbine power (MW)
compCl = [0.11 0.43];       % lift coefficients: [b0(-); b1(deg^-1)]
compCd = [0.0282 0.0042...
    5.3613e-04];  % drag coeff.: [a0(-); a1(deg^-1); a2(deg^-2)]
initTAS = 131.9263;         % initial true airspeed (m/s)
g = 9.81;                   % acceleration due to gravity (m/s^2)
S = 77.3;                   % characteristic surface area (m^2)

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
N = length(t)+1 - timestep;
Phi = ones(N, 1);

% Get drive power coeff. in P_drv = eta(:,1).*m.^2 + eta(:,2).*m + eta(:,3)
[eta2, eta1, eta0] = f_drv(v, v_dot, gamma, gamma_dot, compCd, compCl,...
    g, S, rho_air); % (MW per shaft)
eta = [eta2(timestep:end), eta1(timestep:end), eta0(timestep:end)];

% gt fuel coeff. in ph = beta1*Peng + beta0 (1st order fit)
beta = [0.0821 0.0327].*ones(N, 2); %  (beta = [beta1; beta0])

% e.m. loss map coeff. in Pc = alpha2*P_gt^2 + alpha1*P_gt + alpha0
alpha = [0 1 0].*ones(N, 3); %  (alpha = [alpha2; alpha1; alpha0])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Power constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Bounds on variables
% Bounds on electric motor input (electrical) power P_c
P_c_up = h(P_em_up, alpha);

% Bounds on fuel consumption
phi_low = f(P_gt_low*Phi, beta);
phi_up = f(P_gt_up*Phi, beta);

% Estimate bounds for P_b
P_b_up_est1 = inv_f_P_b(P_c_up, U_max, R_max);
P_b_up_est2 = U_max^2/(2*R_max)*10^-6;
P_b_up_est = min(P_b_up_est1, P_b_up_est2*Phi);

%% scale factors for parameters
E_scale = 1.4875e3;
U_scale = 800;
P_b_scale = 3.6328125;
phi_scale = 0.4432;
M_scale = 42000;
R_scale = 0.0350;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Battery Coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

voc = readmatrix('data/Voc_normalised.csv');
R = readmatrix('data/Rbat_normalised.csv');

if open_circuit_deg == 0
    v_coeff = 1;
else
    v_coeff = polyfit(voc(:,1), voc(:,2), open_circuit_deg);
end

if resistance_deg == 0
    r_coeff = 1;
else
    r_coeff = polyfit(R(:,1), R(:,2), open_circuit_deg);
end

clear voc R;

%% Initialise values
E_0 = min(E_init)*Phi/E_scale;
U_0 = U_max*Phi/U_scale;
m_0 = M_init*Phi/M_scale;
R_0 = R_max*Phi/R_scale;

if timestep == 1
    x0 = repmat([m_0; phi_up/phi_scale; E_0; P_b_up_est/P_b_scale; R_0; U_0], K, 1);
end

%% we take the optimisation variable as x = [m, phi, E, P_b, R, U]*K;
% The callback functions.
funcs.objective         = @objective;
funcs.constraints       = @constraintsMPC;
funcs.gradient          = @gradipopt;
funcs.jacobian          = @jacobianMPC;
funcs.jacobianstructure = @jacobiansparseMPC;

% Lower bound on the variables.
options.lb = repmat([-Inf*ones(N,1); phi_low/phi_scale; E_low*ones(N,1)/E_scale; ...
    -Inf*ones(N,1); - Inf*ones(N,1); -Inf*ones(N,1)], K, 1);
% Upper bound on the variables.
options.ub = repmat([Inf*ones(N,1); phi_up/phi_scale; E_up*ones(N,1)/E_scale; ...
    Inf*ones(N,1); Inf*ones(N,1); Inf*ones(N,1)], K, 1);
% Lower bounds on the constraint functions.
options.cl = [repmat([zeros(4*N,1); -Inf*ones(2*N,1); zeros(2*N,1)], K, 1); zeros(K-1,1)];
% Upper bounds on the constraint functions.
options.cu = [repmat([zeros(N,1); Inf*ones(N,1); zeros(5*N,1); Inf*ones(N,1)], K, 1); zeros(K-1,1)];

% Set the IPOPT options.
options.ipopt.mu_strategy = 'adaptive';
options.ipopt.tol         = 1e-10;
options.ipopt.print_level = 3;
options.ipopt.hessian_approximation = 'limited-memory';
options.auxdata = {N, alpha, beta, eta, U_max, R_max, SOC_max, ... %(1-7)
    M_init, E_init, P_em_up, v_coeff, r_coeff, d, n_mot, ...       %(8-14)
    E_scale, U_scale, P_b_scale, phi_scale, M_scale, R_scale, ...  %(15-20)
    P_em_low, K, w_k};                                             %(21-22)
%% Run IPOPT.
[x, info] = ipopt_auxdata(x0,funcs,options);

%% Write data for IPOPT
U = f_ubat(E_est*E_scale/SOC_max, U_max, v_coeff);
R = f_rbat(E_est*E_scale/SOC_max, R_max, r_coeff);
P_drv = f_P_drv(M_init, eta(1,:));
if info.status == 0
    P_b = P_b_scale*x(3*N+1);
    P_c = inv_f_P_b(P_b, U, R)*10^6;
    phi = max(phi_low(1), f_phi(P_drv, P_b, U, R, alpha(1,:), beta(1,:))/(n_mot));
else
    x=x0;
    P_b = P_b_scale*x(3*N+1);
    P_c = inv_f_P_b(P_b, U, R)*10^6;
    phi = max(phi_low(1), f_phi(P_drv, P_b, U, R, alpha(1,:), beta(1,:))/(n_mot));
end

end

function o = objective(x, auxdata)
N = cell2mat(auxdata(1));
K = cell2mat(auxdata(22));
o = -mean(x(N:N*6:6*(K-1)*N+N));
end

function g = gradipopt(~, auxdata)
N = cell2mat(auxdata(1));
K = cell2mat(auxdata(22));
g = 1/K*repmat([zeros(N-1,1); -1; zeros(N*5,1)], K, 1);
end