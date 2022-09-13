function C = calcC(SOC, P_c, U_coeff, R_coeff)
% Function: Compute partial derivative of the measurement function h with
% respect to SOC
% Input: 
%        - SOC: state of charge (dimensionless)
%        - P_c: power demand of electric motor (W)
%        - SOC_max: maximum state of charge (J)
%        - U_coeff/R_coeff: polynomial coefficients for open-circuit
%        voltage, battery internal resistance
% Output: partial h/ partial SOC

U = polyval(U_coeff, SOC);
dU = polyval(polyder(U_coeff),SOC);
R = polyval(R_coeff, SOC);
dR = polyval(polyder(R_coeff), SOC);
X = sqrt(U^2 - 4*R*P_c);

C = [dR*(-1/(2*R^2)*(U-X) + P_c/(R*X)) + dU/(2*R)*(1-U/X);
    dR*(2*P_c/(U-X)-4*R*P_c^2/(X*(U-X)^2)) - dU*(2*R*P_c)*(1-U/X)/(U-X)^2];
end