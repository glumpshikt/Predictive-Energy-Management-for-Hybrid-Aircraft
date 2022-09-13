%% constraint function for IPOPT (MPC)

function c = constraintsMPC(x, auxdata)
N = cell2mat(auxdata(1));
alpha = cell2mat(auxdata(2));
beta= cell2mat(auxdata(3));
eta = cell2mat(auxdata(4));
U_max = cell2mat(auxdata(5));
R_max = cell2mat(auxdata(6));
SOC_max = cell2mat(auxdata(7));
M_init = cell2mat(auxdata(8));
E_init = cell2mat(auxdata(9));
P_em_up = cell2mat(auxdata(10));
v_coeff = cell2mat(auxdata(11));
r_coeff = cell2mat(auxdata(12));
d = cell2mat(auxdata(13));
n_mot = cell2mat(auxdata(14));
E_scale = cell2mat(auxdata(15));
U_scale = cell2mat(auxdata(16));
P_b_scale = cell2mat(auxdata(17));
phi_scale = cell2mat(auxdata(18));
M_scale = cell2mat(auxdata(19));
R_scale = cell2mat(auxdata(20));
P_em_low = cell2mat(auxdata(21));
K = cell2mat(auxdata(22));
w_k = cell2mat(auxdata(23));

c = zeros(8*N*K+K-1,1);
for k=1:K
    m = x(1+(k-1)*6*N: N+(k-1)*6*N);
    phi = x(N+1+(k-1)*6*N: 2*N+(k-1)*6*N);
    E = x(2*N+1+(k-1)*6*N: 3*N+(k-1)*6*N);
    P_b = x(3*N+1+(k-1)*6*N: 4*N+(k-1)*6*N);
    R = x(4*N+1+(k-1)*6*N: 5*N+(k-1)*6*N);
    U = x(5*N+1+(k-1)*6*N: 6*N+(k-1)*6*N);

    P_drv = f_P_drv(m*M_scale, eta);

    %constraint 1: mass <-> fuel equality
    c(1+(k-1)*8*N: N+(k-1)*8*N) = d*n_mot*cumsum([0; phi(1:end-1)])*phi_scale/M_scale + m - M_init*ones(N,1)/M_scale;
    %constraint 2: c(x) = phi - f_phi(m, P_b, U, R)
    c(N+1+(k-1)*8*N: 2*N+(k-1)*8*N) ...
        = phi - f_phi(P_drv, P_b*P_b_scale, U*U_scale, R*R_scale, alpha, beta)/(phi_scale * n_mot);
    %constraint 3: SOC <-> P_b equality
    c(2*N+1+(k-1)*8*N: 3*N+(k-1)*8*N) ...
        = d*cumsum([0; P_b(1:end-1)])*P_b_scale/E_scale + E + cumsum(w_k(k,:)).'/E_scale - E_init(k)*ones(N,1)/E_scale;
    %constraint 4: SOC <-> V_oc
    c(3*N+1+(k-1)*8*N: 4*N+(k-1)*8*N) ...
        = U - f_ubat(E*E_scale/SOC_max, U_max, v_coeff)/U_scale;
    %constraint 5: P_b upperlim 1
    c(4*N+1+(k-1)*8*N: 5*N+(k-1)*8*N) ...
        = f_P_b_lim1(P_b*P_b_scale, U*U_scale, R*R_scale, alpha, P_em_up);
    %constraint 6: P_b upperlim 2
    c(5*N+1+(k-1)*8*N: 6*N+(k-1)*8*N) ...
        = f_P_b_upperlim2(P_b*P_b_scale, U*U_scale, R*R_scale);
    %constraint 7: SOC <-> Resistance
    c(6*N+1+(k-1)*8*N: 7*N+(k-1)*8*N) ...
        = R - f_rbat(E*E_scale/SOC_max, R_max, r_coeff)/R_scale;
    %constraint 8: P_b lowerlim
    c(7*N+1+(k-1)*8*N: 8*N+(k-1)*8*N) ...
        = f_P_b_lim1(P_b*P_b_scale, U*U_scale, R*R_scale, alpha, P_em_low);
end

%constraint 9: P_b_1,1 = P_b_2,1 = ... = P_b,k,1
% for k=1:K-1
%     c(8*K*N+k) = x(3*N+1) - x(3*N+1+k*6*N);
% end
end