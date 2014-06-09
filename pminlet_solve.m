function res = pminlet_solve(Min,Mout,n_div,H0,phi0,dtheta)
					  
	M = linspace(Min,Mout,n_div)';
	am =  atan(1./sqrt(M.^2-1));

	theta(1,1) = dtheta;
	r(1,1) = 1;
	phi(1,1) = deg2rad(phi0);
	for i = 2:n_div-1
		ds_dM(i,1) =  (-am(i+1)+2*am(i)-am(i-1))/(M(i)-M(i-1));
		theta(i,1) = theta(1,1);
		dr(i,1) = (sin(am(i-1,1)+theta(i-1,1))/sin(am(i,1))-1)*r(i-1,1);
		r(i,1) = r(i-1,1)+dr(i,1);
		dphi(i,1) = am(i,1)-(am(i-1,1)+theta(i-1,1));
		phi(i,1) = phi(i-1,1)-dphi(i,1);
	end
	x = r.*cos(phi);
	y = r.*sin(phi);
	
	for i = 2:n_div-1
		H(i,1) = atan((y(i,1)-y(i-1,1))/(x(i,1)-x(i-1,1)));
	end
	res = H(2,1) - deg2rad(H0);
end