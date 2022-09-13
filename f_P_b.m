function P_b = f_P_b(P_c, U, R)
% Function: Compute battery power from input power to electric motor
% Input:
%        - P_c: Electric motor input power (MW)
%        - U: Battery open-circuit voltage (V)
%        - R: Battery internal resistance (Ohms)
% Output: P_b: Power drawn from battery (MW)

P_b= U.^2./(2*R).*(1-sqrt(1-4*10^6*R.*P_c./(U.^2))); 
end

