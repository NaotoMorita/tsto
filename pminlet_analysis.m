for iter_m = 1
	for iter_aoa = 1
	Mmain_off = 6
	aoa_body = 2
	10*(iter_m-1)+iter_aoa

kappa=1.4

%機首予圧縮による圧縮
if -aoa+aoa_body >=0
	beta_off=fsolve(@(beta)tbm_function(deg2rad(-aoa+aoa_body),Mmain_off,beta,1.4),0.0001);
	[Min_off p02_p01_off p2_p1_off  a2_a1_off T2_T1_off rho2_rho1_off]=pmvalue(Mmain_off,beta_off,1.4)
else
	Min_off = fsolve(@(M2)pmsolve(Mmain_off,M2,-deg2rad(-aoa+aoa_body),kappa),Mmain_off)
	rho2_rho1_off = ((1+(kappa-1)/2* Mmain_off^2)/(1+(kappa-1)/2*Min_off^2))^(1/(kappa-1))
	a2_a1_off = ((1+(kappa-1)/2* Mmain_off^2)/(1+(kappa-1)/2*Min_off^2))^(1/2)
	p02_p01_off = 1;
	p2_p1_off =  ((1+(kappa-1)/2* Mmain_off^2)/(1+(kappa-1)/2*Min_off^2))^(kappa/(kappa-1))
	T2_T1_off =  ((1+(kappa-1)/2* Mmain_off^2)/(1+(kappa-1)/2*Min_off^2))
end
%第1ランプ
if -H0+aoa != 0
	beta_1ramp_off=fsolve(@(beta)tbm_function(deg2rad(-H0+aoa),Min_off,beta,1.4),0.0001);
	[Min_1ramp_off p02_p01_1ramp_off p2_p1_1ramp_off a2_a1_1ramp_off T2_T1_1ramp_off rho2_rho1_1ramp_off]=pmvalue(Min_off,beta_1ramp_off,1.4)
else
	beta_1ramp_off = pi/2;
	Min_1ramp_off = Min_off;
	p02_p01_1ramp_off = 1;
	p2_p1_1ramp_off = 1;
	a2_a1_1ramp_off = 1;
	T2_T1_1ramp_off =1;
	rho2_rho1_1ramp_off=1;

end

%等エントロピー圧縮
M_off(1,1) = Min_1ramp_off; 
p_off(1,1) = p0*p2_p1_1ramp_off*p2_p1_off;
for i = 2: rows(theta)
	M_off(i,1) = fsolve(@(M2)pmsolve(M_off(i-1,1),M2,theta(i),1.4),M_off(i-1,1));
	dp_off(i,1) = -theta(i,1)*1.4*M_off(i,1)^2/sqrt(M_off(i,1)^2-1)*p_off(i-1,1);
	p_off(i,1) = p_off(i-1,1)+dp_off(i,1);
end
rho2_rho1_S_off = ((1+(kappa-1)/2* Min_1ramp_off^2)/(1+(kappa-1)/2*M_off(rows(M_off))^2))^(1/(kappa-1))
a2_a1_S_off = ((1+(kappa-1)/2* Min_1ramp_off^2)/(1+(kappa-1)/2*M_off(rows(M_off))^2))^(1/2)

beta_Tshock_off = fsolve(@(beta)tbm_function(-atan(tan_slast),M_off(rows(M_off),1),beta,1.4),0.001);
shock_angle = beta_Tshock_off+H(rows(H),1)
[Mout_throat_off p02_p01_throat_off p2_p1_throat_off  a2_a1_throat_off T2_T1_throat_off rho2_rho1_throat_off] = pmvalue(M_off(rows(M_off),1),beta_Tshock_off,1.4);
TPR_off =  p02_p01_off*p02_p01_throat_off
MCR_off = Mout_throat_off/Min_off*a2_a1_off*a2_a1_1ramp_off*a2_a1_throat_off*a2_a1_S_off*rho2_rho1_off*rho2_rho1_1ramp_off*rho2_rho1_throat_off*rho2_rho1_S_off/cos(deg2rad(-H0))*ys/y_nose
p_pin_off = p2_p1_throat_off*p_off(rows(p_off),1)/p0

%[F_off Isp_off] = scram_thrust(MCR_off,TPR_off,p_pin_off,Tmc,Mmain_off,216,5531,0.088,max(y_nose)/1000*2.02,ys/1000*2.02,3.66*2.02,1)

M_anly(iter_m) = Mmain_off
aoa_anly(iter_aoa) = aoa_body
TPR_anly(iter_m,iter_aoa) = TPR_off
MCR_anly(iter_m,iter_aoa) = MCR_off
p_pin_anly(iter_m,iter_aoa) = p_pin_off
F_anly(iter_m,iter_aoa) =F_off
Isp_anly(iter_m,iter_aoa) =Isp_off
end

end