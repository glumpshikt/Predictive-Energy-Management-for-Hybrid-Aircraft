%% jacobian for IPOPT (MPC)

function j = jacobianMPC(x, auxdata)
N = cell2mat(auxdata(1));
alpha = cell2mat(auxdata(2));
beta= cell2mat(auxdata(3));
eta = cell2mat(auxdata(4));
U_max = cell2mat(auxdata(5));
R_max = cell2mat(auxdata(6));
SOC_max = cell2mat(auxdata(7));
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
K = cell2mat(auxdata(22));

m = x(1:N);

I = speye(N);
Psi = sparse(tril(ones(N), -1));

j = zeros(8*N*K+K-1, 6*N*K);

for k=1:K
    m = x(1+(k-1)*6*N: N+(k-1)*6*N);
    E = x(2*N+1+(k-1)*6*N: 3*N+(k-1)*6*N);
    P_b = x(3*N+1+(k-1)*6*N: 4*N+(k-1)*6*N);
    R = x(4*N+1+(k-1)*6*N: 5*N+(k-1)*6*N);
    U = x(5*N+1+(k-1)*6*N: 6*N+(k-1)*6*N);
    
    %constraint 1: mass <-> fuel equality
    j(1+(k-1)*8*N: N+(k-1)*8*N, 1+(k-1)*6*N:2*N+(k-1)*6*N) = [I Psi*d*n_mot*phi_scale/M_scale];
    
    %constraint 2: c(x) = phi - f_phi(m, P_b, U)
    c2dU = beta(:,1)*U_scale/(phi_scale * n_mot) ...
        .* dinv_h(inv_f_P_b(P_b*P_b_scale, U*U_scale, R*R_scale),alpha) ...
        .* dinv_f_P_bU(P_b*P_b_scale, R*R_scale, U*U_scale);
    c2dP_b = beta(:,1)*P_b_scale/(phi_scale * n_mot) ...
        .* dinv_h(inv_f_P_b(P_b*P_b_scale, U*U_scale, R*R_scale),alpha) ...
        .* dinv_f_P_bP_b(P_b*P_b_scale, R*R_scale, U*U_scale);
    c2dphi = ones(N,1);
    c2dm = (-1)* beta(:,1)*M_scale/(phi_scale * n_mot) ...
        .* df_P_drv(m*M_scale, eta);
    c2dR = beta(:,1)*R_scale/(phi_scale * n_mot) ...
        .* dinv_h(inv_f_P_b(P_b*P_b_scale, U*U_scale, R*R_scale),alpha) ...
        .* dinv_f_P_bR(P_b*P_b_scale, R*R_scale, U*U_scale);
    j(N+1+(k-1)*8*N: 2*N+(k-1)*8*N, 1+(k-1)*6*N:2*N+(k-1)*6*N) ...
        = [spdiags(c2dm,0,N,N) spdiags(c2dphi,0,N,N)];
    j(N+1+(k-1)*8*N: 2*N+(k-1)*8*N, 3*N+1+(k-1)*6*N: 6*N+(k-1)*6*N) ...
        = [spdiags(c2dP_b,0,N,N) spdiags(c2dR,0,N,N) spdiags(c2dU,0,N,N)];
    clear c2dU c2dP_b c2dphi c2dm c2dR;

    %constraint 3: SOC <-> P_b equality
    j(2*N+1+(k-1)*8*N: 3*N+(k-1)*8*N, 2*N+1+(k-1)*6*N: 4*N+(k-1)*6*N) ...
        = [I Psi*d*P_b_scale/E_scale];

    %constraint 4: SOC <-> V_oc
    c4dE = -df_ubat(E*E_scale/SOC_max, U_max, v_coeff)*E_scale/(SOC_max*U_scale);
    j(3*N+1+(k-1)*8*N: 4*N+(k-1)*8*N, 2*N+1+(k-1)*6*N: 3*N+(k-1)*6*N) ...
        = spdiags(c4dE,0,N,N);
    j(3*N+1+(k-1)*8*N: 4*N+(k-1)*8*N, 5*N+1+(k-1)*6*N: 6*N+(k-1)*6*N) ...
        = I;
    clear c4dE;

    %constraint 5: P_b upperlim 1
    c5dU = U_scale*dinv_h(inv_f_P_b(P_b*P_b_scale, U*U_scale, R*R_scale),alpha) ...
        .* dinv_f_P_bU(P_b*P_b_scale, R*R_scale, U*U_scale);

    c5dP_b = P_b_scale*dinv_h(inv_f_P_b(P_b*P_b_scale, U*U_scale, R*R_scale),alpha) ...
        .* dinv_f_P_bP_b(P_b*P_b_scale, R*R_scale, U*U_scale);

    c5dR = R_scale*dinv_h(inv_f_P_b(P_b*P_b_scale, U*U_scale, R*R_scale),alpha) ...
        .* dinv_f_P_bR(P_b*P_b_scale, R*R_scale, U*U_scale);
    j(4*N+1+(k-1)*8*N: 5*N+(k-1)*8*N, 3*N+1+(k-1)*6*N: 6*N+(k-1)*6*N) ...
        = [spdiags(c5dP_b,0,N,N) spdiags(c5dR,0,N,N) spdiags(c5dU,0,N,N)];
    clear c5dU c5dP_b c5dR;

    %constraint 6: P_b upperlim 2
    c6dR = 0.5*U.^2./(R.^2)*10^-6*R_scale;
    c6dU = -U./R*10^-6*U_scale;
    j(5*N+1+(k-1)*8*N: 6*N+(k-1)*8*N, 3*N+1+(k-1)*6*N: 6*N+(k-1)*6*N) ...
        = [P_b_scale*I spdiags(c6dR,0,N,N) spdiags(c6dU,0,N,N)];
    clear c6dR c6dU;

    %constraint 7: SOC <-> R
    c7dE = -df_rbat(E*E_scale/SOC_max, R_max, r_coeff)*E_scale/(SOC_max*R_scale);
    j(6*N+1+(k-1)*8*N: 7*N+(k-1)*8*N, 2*N+1+(k-1)*6*N: 3*N+(k-1)*6*N) ...
        = spdiags(c7dE,0,N,N);
    j(6*N+1+(k-1)*8*N: 7*N+(k-1)*8*N, 4*N+1+(k-1)*6*N: 5*N+(k-1)*6*N) ...
        = I;
    clear c7dE;
    
    %constraint 8: P_b lowerlim
    j(7*N+1+(k-1)*8*N: 8*N+(k-1)*8*N, 3*N+1+(k-1)*6*N: 6*N+(k-1)*6*N) = ...
        j(4*N+1+(k-1)*8*N: 5*N+(k-1)*8*N, 3*N+1+(k-1)*6*N: 6*N+(k-1)*6*N);
end

for k=1:K-1
    j(8*K*N+k, 3*N+1) = 1;
    j(8*K*N+k, 3*N+1+k*6*N) = -1;
end

j = sparse(j);
end