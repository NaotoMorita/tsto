function res = tbm_function(theta, mach, beta, kappa)
	res = tan(theta) - ((mach ^ 2 * sin(beta) ^ 2 - 1) * cot(beta) / (1 + (1 / 2 * (kappa + 1) - sin(beta) ^ 2) * mach^2));
endfunction