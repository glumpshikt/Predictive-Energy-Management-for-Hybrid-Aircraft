function out = dinv_h(P_c, alpha)
%inverse electric motor map
out = 0;
if (alpha(1, 3) == 0)  % linear map
    out = 1./alpha(:,2);
elseif(alpha(1, 3) ~= 0)  %quadratic map
    out = (-alpha(:,2) + sqrt(alpha(:,2).^2 ...
        -4*alpha(:,3).*(alpha(:,1)-P_c)))./(2*alpha(:,3));
end
end

