function out = dinv_f_P_bU(P_b, R, U)
% Function: Compute partial derivative of dinv_f_P_b with respect to U
% Input:  
%         - P_c: effective electrical power (MW)
%         - R: Battery internal resistance (Ohms)
%         - U: Open circuit voltage (V)
% Output: 
%         - out: dinv_f_P_b/dU

out = 2*R.*(P_b.^2)./(U.^3)*10^6;
end

