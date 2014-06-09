function res = pminlet_k_solve(Min,Mout,n_div,H0,hin,p0,k)
	theta0 = prandtl_meyer_function(Min,1.4)-prandtl_meyer_function(Mout,1.4); %deg
	dtheta = -theta0/(n_div) * k;
	%phi0 = 170
	%res = pminlet_solve(Min,Mout,n_div,H0,phi0)
	%phi0 = fsolve(@(phi0)pminlet_solve(Min,Mout,n_div,H0,phi0,deg2rad(dtheta)),179)

	M = linspace(Min,Mout,n_div)';
	am =  atan(1./sqrt(M.^2-1));
	phi0 =  180 - rad2deg(am(1,1)) + H0;
	theta(1,1) = deg2rad(dtheta);
	r(1,1) = 1;
	p(1,1) = p0;
	phi(1,1) = deg2rad(phi0);
	for i = 2:n_div-1
		theta(i,1) = -prandtl_meyer_function(M(i,1),1.4)+prandtl_meyer_function(M(i-1,1),1.4)
		dp(i,1) = -theta(i,1)*1.4*M(i,1)^2/sqrt(M(i,1)^2-1)*p(i-1,1);
		p(i,1) = p(i-1,1)+dp(i,1);
		dr(i,1) = (sin(am(i-1,1)+theta(i-1,1))/sin(am(i,1))-1)*r(i-1,1);
		r(i,1) = r(i-1,1)+dr(i,1);
		dphi(i,1) = am(i,1)-(am(i-1,1)+theta(i-1,1));
		phi(i,1) = phi(i-1,1)-dphi(i,1);
	end
	H(1,1) = deg2rad(H0);
	x = r.*cos(phi);
	y = r.*sin(phi);

	for i = 2:n_div-1
		H(i,1) = atan((y(i)-y(i-1))/(x(i)-x(i-1)));
	end

	%インテークの相似拡大
	x_intake = hin./max(y) .*x;
	y_intake = hin./max(y) .*y;
	r0 = sqrt(x_intake(1,1)^2+y_intake(1,1)^2);

	%接続部は5次関数

	%衝撃波角計算
	tan_slast = (y_intake(rows(y),1)- y_intake(rows(y)-1,1))/(x_intake(rows(x),1)- x_intake(rows(x)-1,1));
	beta = fsolve(@(beta)tbm_function(-atan(tan_slast),M(rows(M),1),beta,1.4),0.001);
	shock_angle = beta+H(rows(H),1);
	[mach2 p02_p01 p2_p1] = pmvalue(M(rows(M),1),beta,1.4);
	MCR = 1;
	TPR = p02_p01;
	mach2;
	p_pin = p2_p1*p(rows(p),1)/p0;
	xo = x_intake(rows(x_intake),1);
	yo = y_intake(rows(x_intake),1);

	x_sholder = (-xo*tan_slast+yo)/(-tan_slast+tan(shock_angle));
	xs = x_sholder;
	ys = tan(shock_angle)*x_sholder;
	ht = ys;
	x_sp = linspace(xo+0.00001,xs,500)';
	y_sp = tan_slast.*(x_sp-xo)+yo;

	tan_ain = -1/tan_slast;
	x_in = (-xo*tan_slast+yo)/(-tan_slast+tan_ain);
	y_in = tan_ain*x_in;
	Ain = sqrt(x_in^2+y_in^2);
	At_Ain = ys/Ain;
	beta =rad2deg(beta);
	theta0;
	shock_theta =- rad2deg(tan_slast);

	res = shock_theta-theta0+H0
end