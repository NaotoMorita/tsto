function nu = prandtl_meyer_function(M,kappa)
	nu = sqrt((kappa+1)/(kappa-1))*atan(sqrt((kappa-1)/(kappa+1)*(M^2-1)))-atan(sqrt(M^2-1));
	nu =nu*180/pi;
end