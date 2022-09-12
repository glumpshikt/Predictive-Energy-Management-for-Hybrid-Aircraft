function r_bat = df_rbat(SOC, R_max, coeff)
% Function: Compute r_bat map (battery internal resistance map)
% Input:
%        - SOC: State-of-Charge (dimensionless)
%        - R_max: Maximum resistance (Ohms)
%        - coefficients of polynomial
% Output: df_rbat/dSOC

r_bat = R_max*polyval(polyder(coeff), SOC);
end

