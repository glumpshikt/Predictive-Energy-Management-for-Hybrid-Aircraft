function P_b_upperlim2 = f_P_b_upperlim2(P_b, U, R)
% Function: Compute P_b constraint 2 map
% Input:
%        - P_b: Battery chemical power (MW)
%        - U: Battery open-circuit voltage (V)
%        - R: Battery internal resistance (Ohms)
% Output: P_b constraint 2

P_b_upperlim2 = P_b - U.^2./(2*R)*10^-6;

end

