%% Function to update state f
function [SOC_updated, U_updated, R_updated] = stateupdate(SOC, w, ...
    P_c, d, SOC_max, U_coeff, R_coeff)
U = polyval(U_coeff, SOC);
R = polyval(R_coeff, SOC);
SOC_updated = SOC - 1/(2*R)*(U-sqrt(U^2-4*R*P_c))*U*d/(SOC_max)+w;
U_updated = polyval(U_coeff, SOC_updated);
R_updated = polyval(R_coeff, SOC_updated);
end