% Calculate matrix of partial f/ partial SOC

function A = calcA(SOC, P_c, d, SOC_max, U_coeff, R_coeff)
U = polyval(U_coeff, SOC);
dU = polyval(polyder(U_coeff),SOC);
R = polyval(R_coeff, SOC);
dR = polyval(polyder(R_coeff), SOC);
X = sqrt(U^2 - 4*R*P_c);

A = 1+(d/SOC_max)*dR*(U/(2*R^2)*(U-X) - U*P_c/(R*X)) ...
    - (d/SOC_max)*dU*(1/(2*R)*(U-X)+U/(2*R)*(1-U/X));

end