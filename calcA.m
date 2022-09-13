function A = calcA(SOC, P_c, d, SOC_max, U_coeff, R_coeff)
% Function: Compute partial derivative of the state update function f with
% respect to SOC
% Input: 
%        - SOC: state of charge (dimensionless)
%        - P_c: power demand of electric motor (W)
%        - d: duration of timestep (s)
%        - SOC_max: maximum state of charge (J)
%        - U_coeff/R_coeff: polynomial coefficients for open-circuit
%        voltage, battery internal resistance
% Output: partial f/ partial SOC

U = polyval(U_coeff, SOC);
dU = polyval(polyder(U_coeff),SOC);
R = polyval(R_coeff, SOC);
dR = polyval(polyder(R_coeff), SOC);
X = sqrt(U^2 - 4*R*P_c);

A = 1+(d/SOC_max)*dR*(U/(2*R^2)*(U-X) - U*P_c/(R*X)) ...
    - (d/SOC_max)*dU*(1/(2*R)*(U-X)+U/(2*R)*(1-U/X));

end