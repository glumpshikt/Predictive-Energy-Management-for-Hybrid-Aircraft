function D = calcD(SOC, P_c, U_coeff, R_coeff)
% Function: Compute partial derivative of the measurement function h with
% respect to SOC
% Input: 
%        - SOC: state of charge (dimensionless)
%        - P_c: power demand of electric motor (W)
%        - SOC_max: maximum state of charge (J)
%        - U_coeff/R_coeff: polynomial coefficients for open-circuit
%        voltage, battery internal resistance
% Output: partial h/ partial P_c

U = polyval(U_coeff, SOC);
R = polyval(R_coeff, SOC);
X = sqrt(U^2 - 4*R*P_c);

D = [1/X; 2*R/(U-X)-4*R^2*P_c/(X*(U-X)^2)];
end