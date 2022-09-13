function out = df_ubat(SOC, U_max, coeff)
% Function: Compute derivative of u_bat map (Open Circuit Voltage map)
% Input:
%        - SOC: State-of-Charge (dimensionless)
%        - U_max: Maximum open circuit voltage (V)
%        - coeff: coefficients of polynomial
% Output: df_ubat/dSOC: 

out = U_max*polyval(polyder(coeff), SOC);
end

