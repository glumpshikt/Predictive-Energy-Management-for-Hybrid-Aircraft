function r_bat = f_rbat(SOC, R_max, coeff)
% Function: Compute r_bat map (battery internal resistance map)
% Input:
%        - SOC: State-of-Charge (dimensionless)
%        - R_max: Maximum resistance (Ohms)
%        - coefficients of polynomial
% Output: u_bat: open circuit voltage (V)
r_bat = R_max*polyval(coeff, SOC);
end

