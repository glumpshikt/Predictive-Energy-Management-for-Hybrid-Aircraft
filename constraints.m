%% constraint function for IPOPT

function c = constraints(x, auxdata)
N = cell2mat(auxdata(1));
alpha = cell2mat(auxdata(2));
beta= cell2mat(auxdata(3));
eta = cell2mat(auxdata(4));
U_max = cell2mat(auxdata(5));
R_max = cell2mat(auxdata(6));
SOC_max = cell2mat(auxdata(7));
M = cell2mat(auxdata(8));
E_up = cell2mat(auxdata(9));
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

E = x(1:N);
U = x(N+1:2*N);
P_b = x(2*N+1:3*N);
phi = x(3*N+1:4*N);
m = x(4*N+1:5*N);
R = x(5*N+1:6*N);

P_drv = f_P_drv(m*M_scale, eta);

%constraint 1: c(x) = phi - f_phi(m, P_b, U, R)
c = phi - f_phi(P_drv, P_b*P_b_scale, U*U_scale, R*R_scale, alpha, beta)/(phi_scale * n_mot);
%constraint 2: mass <-> fuel equality
c = [c; d*n_mot*cumsum([0; phi(1:end-1)])*phi_scale/M_scale + m - M*ones(N,1)/M_scale];
%constraint 3: SOC <-> P_b equality
c = [c; d*cumsum([0; P_b(1:end-1)])*P_b_scale/E_scale + E - E_up*ones(N,1)/E_scale];
%constraint 4: SOC <-> V_oc
c = [c; U - f_ubat(E*E_scale/SOC_max, U_max, v_coeff)/U_scale];
%constraint 5: P_b upperlim 1
c = [c; f_P_b_lim1(P_b*P_b_scale, U*U_scale, R*R_scale, alpha, P_em_up)];
%constraint 6: P_b upperlim 2
c = [c; f_P_b_upperlim2(P_b*P_b_scale, U*U_scale, R*R_scale)];
%constraint 7: SOC <-> Resistance
c = [c; R - f_rbat(E*E_scale/SOC_max, R_max, r_coeff)/R_scale];
%constraint 8: P_b lowerlim
c = [c; f_P_b_lim1(P_b*P_b_scale, U*U_scale, R*R_scale, alpha, P_em_low)];
end