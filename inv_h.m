function P_em = inv_h(P_c, alpha)
% Function: inverse input/output map for the electric motor.
% Compute the electric motor power P_em (MW) for a input power P_c (MW)
% Input:  
%         - P_c: effective electrical power (MW)
%         - alpha : coefficients of map
% Output: 
%         - P_em: electric motor power (MW)

P_em = 0; 
if (alpha(1, 3) == 0)  % linear map
    P_em = (P_c - alpha(:,1))./alpha(:,2);
elseif(alpha(1, 3) ~= 0)  % quadratic map
    rho = alpha(:, 2).^2 - 4*alpha(:, 3).*(alpha(:, 1) - P_c);
    P_em = (-alpha(:, 2) + sqrt(rho))./(2*alpha(:,3));
end 
end

