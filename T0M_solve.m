function res = T0M_solve(T0,T,M)
	kappa = 1.4;
	res = T*(1+(kappa-1)/2*M^2)-T0;
end