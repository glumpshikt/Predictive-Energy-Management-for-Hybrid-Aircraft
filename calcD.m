% Calculate matrix of partial h/ partial SOC


function D = calcD(SOC, P_c, U_coeff, R_coeff)

U = polyval(U_coeff, SOC);
R = polyval(R_coeff, SOC);
X = sqrt(U^2 - 4*R*P_c);

D = [1/X; 2*R/(U-X)-4*R^2*P_c/(X*(U-X)^2)];
end