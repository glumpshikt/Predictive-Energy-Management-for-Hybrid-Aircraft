function P_b = f_P_b(P_c, U, R)

P_b= U.^2./(2*R).*(1-sqrt(1-4*10^6*R.*P_c./(U.^2))); 
end

