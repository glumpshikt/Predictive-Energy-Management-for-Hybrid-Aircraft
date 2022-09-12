function P_drv = f_P_drv(m, eta)
% Function: Compute P_drv map (Total Power Requirement map)
% Input: 
%        - m: mass (kg)
%        - m_drv, n_drv, q_drv (map coefficients)
% Output: P_drv: Total Power (MW)

P_drv = eta(:,1).*(m.^2) + eta(:,2).*m + eta(:,3); 
end

