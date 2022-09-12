function F_phi = f_phi(P_drv, P_b, U, R, alpha, beta)
% Function: Compute f_phi map (rate of change of fuel map)
% Input: 
%        - P_drv: drive power per system (MW)
%        - P_b: battery chemical power per system (MW)
%        - alpha: em loss map coefficients
%        - beta: fuel map coefficients
%        - R_U: characteristic ratio R_bat/U_bat^2*10^6
% Output: F_phi: rate of change of fuel mass per engine(kg/s)

F_phi = f(P_drv - inv_h(inv_f_P_b(P_b, U, R), alpha), beta); 
end

