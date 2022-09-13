function [I_out, V_out] = statemeasure(SOC, eps1, eps2, P_c, U_coeff, R_coeff)
% Function: Compute measurement at i+1 given state at timestep i
% Input: 
%        - SOC: state of charge (dimensionless)
%        - eps1: Current measurement noise(set to 0 in SOC estimation)
%        - eps2: Voltage measurement noise(set to 0 in SOC estimation)
%        - P_c: power demand of electric motor (W)
%        - SOC_max: maximum state of charge (J)
%        - U_coeff/R_coeff: polynomial coefficients for open-circuit
%        voltage, battery internal resistance
% Output: values of current and voltage
U = polyval(U_coeff, SOC);
R = polyval(R_coeff, SOC);
I_out = 1/(2*R)*(U-sqrt(U^2-4*R*P_c))+eps1;
V_out = 2*R*P_c/(U-sqrt(U^2-4*R*P_c))+eps2;
end