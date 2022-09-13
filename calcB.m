function B = calcB(SOC, P_c, d, SOC_max, U_coeff, R_coeff)
% Function: Compute partial derivative of the state update function f with
% respect to P_c
% Input: 
%        - SOC: state of charge (dimensionless)
%        - P_c: power demand of electric motor (W)
%        - d: duration of timestep (s)
%        - SOC_max: maximum state of charge (J)
%        - U_coeff/R_coeff: polynomial coefficients for open-circuit
%        voltage, battery internal resistance
% Output: partial f/ partial P_c

U = polyval(U_coeff, SOC);
R = polyval(R_coeff, SOC);
X = sqrt(U^2 - 4*R*P_c);

B = -d*U/(SOC_max*X);
end