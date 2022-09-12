function u_bat = f_ubat(SOC, U_max, coeff)
% Function: Compute u_bat map (Open Circuit Voltage map)
% Input:
%        - SOC: State-of-Charge (dimensionless)
%        - U_max: Maximum open circuit voltage (V)
%        - coefficients of polynomial
% Output: u_bat: open circuit voltage (V)
u_bat = U_max*polyval(coeff, SOC);
end

