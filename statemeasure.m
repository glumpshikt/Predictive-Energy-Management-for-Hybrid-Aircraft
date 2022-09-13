%% state measurement function h
function [I_out, V_out] = statemeasure(SOC, eps1, eps2, P_c, U_coeff, R_coeff)
U = polyval(U_coeff, SOC);
R = polyval(R_coeff, SOC);
I_out = 1/(2*R)*(U-sqrt(U^2-4*R*P_c))+eps1;
V_out = 2*R*P_c/(U-sqrt(U^2-4*R*P_c))+eps2;
end