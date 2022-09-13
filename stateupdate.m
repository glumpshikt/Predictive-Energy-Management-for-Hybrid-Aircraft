function [SOC_updated, U_updated, R_updated] = stateupdate(SOC, w, ...
    P_c, d, SOC_max, U_coeff, R_coeff)
% Function: Compute state at timestep i+1 given state at timestep i
% Input: 
%        - SOC: state of charge (dimensionless)
%        - w: state noise (set to 0 in SOC estimation)
%        - P_c: power demand of electric motor (W)
%        - d: duration of timestep (s)
%        - SOC_max: maximum state of charge (J)
%        - U_coeff/R_coeff: polynomial coefficients for open-circuit
%        voltage, battery internal resistance
% Output: values of SOC, U and R in the next timestep

U = polyval(U_coeff, SOC);
R = polyval(R_coeff, SOC);
SOC_updated = SOC - 1/(2*R)*(U-sqrt(U^2-4*R*P_c))*U*d/(SOC_max)+w;
U_updated = polyval(U_coeff, SOC_updated);
R_updated = polyval(R_coeff, SOC_updated);
end