%% jacobian for IPOPT

function j = jacobian(x, auxdata)
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

E = x(1:N);
U = x(N+1:2*N);
P_b = x(2*N+1:3*N);
m = x(4*N+1: 5*N);
R = x(5*N+1: 6*N);
I = speye(N);
Zero = sparse(N,N);
Psi = sparse(tril(ones(N), -1));

%constraint 1: c(x) = phi - f_phi(m, P_b, U)
c1dU = beta(:,1)*U_scale/(phi_scale * n_mot) ...
    .* dinv_h(inv_f_P_b(P_b*P_b_scale, U*U_scale, R*R_scale),alpha) ...
    .* dinv_f_P_bU(P_b*P_b_scale, R*R_scale, U*U_scale);
c1dP_b = beta(:,1)*P_b_scale/(phi_scale * n_mot) ...
    .* dinv_h(inv_f_P_b(P_b*P_b_scale, U*U_scale, R*R_scale),alpha) ...
    .* dinv_f_P_bP_b(P_b*P_b_scale, R*R_scale, U*U_scale);
c1dphi = ones(N,1);
c1dm = (-1)* beta(:,1)*M_scale/(phi_scale * n_mot) ...
    .* df_P_drv(m*M_scale, eta);
c1dR = beta(:,1)*R_scale/(phi_scale * n_mot) ...
    .* dinv_h(inv_f_P_b(P_b*P_b_scale, U*U_scale, R*R_scale),alpha) ...
    .* dinv_f_P_bR(P_b*P_b_scale, R*R_scale, U*U_scale);
gc1 = [Zero spdiags(c1dU,0,N,N) spdiags(c1dP_b,0,N,N) ...
    spdiags(c1dphi,0,N,N) spdiags(c1dm,0,N,N) spdiags(c1dR,0,N,N)];
clear c1dU c1dP_b c1dphi c1dm c1dR;

%constraint 2: mass <-> fuel equality
gc2 = [Zero Zero Zero Psi*d*n_mot*phi_scale/M_scale I Zero];

%constraint 3: SOC <-> P_b equality
gc3 = [I Zero Psi*d*P_b_scale/E_scale Zero Zero Zero];

%constraint 4: SOC <-> V_oc
c4dE = -df_ubat(E*E_scale/SOC_max, U_max, v_coeff)*E_scale/(SOC_max*U_scale);
gc4 = [spdiags(c4dE,0,N,N) I Zero Zero Zero Zero];
clear c4dE;

%constraint 5: P_b upperlim 1
c5dU = U_scale*dinv_h(inv_f_P_b(P_b*P_b_scale, U*U_scale, R*R_scale),alpha) ...
    .* dinv_f_P_bU(P_b*P_b_scale, R*R_scale, U*U_scale);

c5dP_b = P_b_scale*dinv_h(inv_f_P_b(P_b*P_b_scale, U*U_scale, R*R_scale),alpha) ...
    .* dinv_f_P_bP_b(P_b*P_b_scale, R*R_scale, U*U_scale);

c5dR = R_scale*dinv_h(inv_f_P_b(P_b*P_b_scale, U*U_scale, R*R_scale),alpha) ...
    .* dinv_f_P_bR(P_b*P_b_scale, R*R_scale, U*U_scale);
gc5 = [Zero spdiags(c5dU,0,N,N) spdiags(c5dP_b,0,N,N) Zero Zero spdiags(c5dR,0,N,N)];
clear c5dU c5dP_b c5dR;

%constraint 6: P_b upperlim 2
c6dR = 0.5*U.^2./(R.^2)*10^-6;
gc6 = [Zero spdiags(-U./R*10^-6,0,N,N)*U_scale I*P_b_scale Zero Zero spdiags(c6dR,0,N,N)*R_scale];
clear c6dR;

%constraint 7: SOC <-> R
c7dE = -df_rbat(E*E_scale/SOC_max, R_max, r_coeff)*E_scale/(SOC_max*R_scale);
gc7 = [spdiags(c7dE,0,N,N) Zero Zero Zero Zero I];
clear c7dE;

j = [gc1; gc2; gc3; gc4; gc5; gc6; gc7; gc5];
end