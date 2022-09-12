function out = dinv_f_P_bP_b(P_b, R, U)
% Input:  
%         - P_c: effective electrical power (MW)
%         - R: Battery internal resistance (Ohms)
%         - U: Open circuit voltage (V)
% Output: 
%         - out: dinv_f_P_b/dR

out = ones(size(P_b)) - 2*R.*P_b./(U.^2)*10^6;
end

