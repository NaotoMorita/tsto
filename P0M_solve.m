function res = P0M_solve(P0,p,M)
	kappa = 1.4;
	res = p*(1+(kappa-1)/2*M^2)^((kappa)/(kappa-1))-P0;
end