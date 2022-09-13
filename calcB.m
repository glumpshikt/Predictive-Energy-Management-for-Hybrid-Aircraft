% Calculate matrix of partial f/ partial P_c

function B = calcB(SOC, P_c, d, SOC_max, U_coeff, R_coeff)

U = polyval(U_coeff, SOC);
R = polyval(R_coeff, SOC);
X = sqrt(U^2 - 4*R*P_c);

B = -d*U/(SOC_max*X);
end