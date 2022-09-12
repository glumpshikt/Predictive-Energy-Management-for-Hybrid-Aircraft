function out = df_P_drv(m, eta)
% Function: Compute P_drv map (Total Power Requirement map)
% Input: 
%        - m: mass (kg)
%        - eta (map coefficients)
% Output: dP_drv/dm

out = 2*eta(:,1).*m + eta(:,2); 
end

