function [mach2 p02_p01 p2_p1 a2_a1 T2_T1 rho2_rho1] = pmvalue(mach,beta,kappa)
	ms = mach * sin(beta);
	mach2     = ((1.0 + (kappa - 1.0) * ms ^2 + mach ^ 4.0 * (1.0 / 4.0 * (kappa + 1.0) ^ 2.0 - kappa * sin(beta) ^ 2.0) * sin(beta) ^ 2.0) / ((1.0 + 1.0 / 2.0 * (kappa - 1.0) * ms ^ 2.0) * (kappa * ms ^ 2.0 - 1.0 / 2.0 * (kappa - 1.0)))) ^ (1.0 / 2.0);
	p2_p1     = (kappa * ms ^ 2.0  - 1.0 / 2.0 * (kappa - 1.0)) / (1.0 / 2.0 * (kappa + 1.0));
	calcbuff = (1 + 2 * kappa / (kappa + 1) * (ms ^ 2 - 1)) ^ (-1 / (kappa-1));
	T2_T1     = (1.0 + 1.0 / 2.0 * (kappa - 1.0) * ms ^ 2.0) * (kappa * ms ^ 2.0 - 1.0 / 2.0 * (kappa - 1.0)) / ( 1.0 / 4.0 * (kappa + 1.0) ^ 2.0 * ms ^ 2.0)
	a2_a1     = T2_T1 ^ (1.0 / 2.0)
	rho2_rho1 = 1.0 / 2.0 * (kappa + 1.0) * ms ^ 2.0 / (1.0 + 1.0 / 2.0 * (kappa - 1.0) * ms ^ 2.0)
        p02_p01 = calcbuff * ((kappa+1)*ms^2/((kappa-1) * ms ^ 2+2)) ^ (kappa/(kappa-1));
end