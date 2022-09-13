% Calculate matrix of partial h/ partial SOC

function C = calcC(SOC, P_c, U_coeff, R_coeff)

U = polyval(U_coeff, SOC);
dU = polyval(polyder(U_coeff),SOC);
R = polyval(R_coeff, SOC);
dR = polyval(polyder(R_coeff), SOC);
X = sqrt(U^2 - 4*R*P_c);

C = [dR*(-1/(2*R^2)*(U-X) + P_c/(R*X)) + dU/(2*R)*(1-U/X);
    dR*(2*P_c/(U-X)-4*R*P_c^2/(X*(U-X)^2)) - dU*(2*R*P_c)*(1-U/X)/(U-X)^2];
end