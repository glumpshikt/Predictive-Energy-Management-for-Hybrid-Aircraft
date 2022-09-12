function P_c = inv_f_P_b(P_b, U, R)
% Function: inverse input/output map for the battery electric bus.
% Compute the effective electrical power P_c (MW) delivered to the em for a
% given input battery chemical power P_b (MW).
% Input:  
%         - P_b: battery chemical power (MW) 
%         - R: Resistance of the battery (Ohms)
%         - U: Voltage of the battery (V)
%         circuit resistance and voltage source respectively 
% Output: 
%         - P_c: effective electrical power (MW)

P_c = P_b - (R./(U.^2).*P_b.^2)*10^6; 
end

