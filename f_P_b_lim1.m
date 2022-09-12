function P_b_lim1 = f_P_b_lim1(P_b, U, R, alpha, P_em_lim)
% Function: Compute P_b constraint 1 map
% Input:
%        - P_b: Battery chemical power (MW)
%        - U: Battery open-circuit voltage (V)
%        - R: Battery internal resistance (Ohms)
% Output: P_b constraint 1 

P_b_lim1 = inv_h(inv_f_P_b(P_b, U, R), alpha) - P_em_lim;

end

