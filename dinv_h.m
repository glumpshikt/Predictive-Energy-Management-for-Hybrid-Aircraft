function P_em = dinv_h(P_c, alpha)
% Function: inverse of Electric motor input-output map
% Input:
%          - P_c:  Electric motor input power (MW)
%          - alpha: vector of map coefficients
% Output:
%          - P_em: Electric motor output power (MW)   
P_em = 0;
if (alpha(1, 3) == 0)           % linear map
    P_em = 1./alpha(:,2);
elseif(alpha(1, 3) ~= 0)        % quadratic map
    P_em = (-alpha(:,2) + sqrt(alpha(:,2).^2 ...
        -4*alpha(:,3).*(alpha(:,1)-P_c)))./(2*alpha(:,3));
end
end

