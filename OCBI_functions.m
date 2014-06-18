1;
function dV =  taylor_maccoll(V,theta,kappa)
	dV(1) = V(2);
	x=(kappa-1)/2* (1-V(1)^2-V(2)^2);
	dV(2) = (V(2)^2*V(1)-2*x*V(1)-x*V(2)/tan(theta)) / (x-V(2)^2);
end

function res = shock_upper(a,x_rec,y_rec)
	res = (1-a)^2-(x_rec^2+(y_rec-a)^2);
end

function dr =  shock_design(r,s,c,s0,R0)
	dr(1) = r(2);
	%r(3)= -c*sin(a*(s0-s))+R0
	dr(2) =1/r(1)*(r(1)^2+2*r(2)^2-((r(1)^2+r(2)^2)^(3/2))*r(3));
	dr(3) = c*cos(1.63*(s0-s));
end

function res = clc_solve(clc,x_rec,y_rec,r_uc,y_uc)
	i=1;
	r(1) = sqrt(y_rec^2+x_rec^2);
	s0lc = atan(y_rec/x_rec);
	slc= linspace(atan(y_rec/x_rec),0,500)';
	rlc(1,:) = [r(1),(y_uc*r(1)*cos(s0lc))/(r(1)-y_uc*sin(s0lc)),r_uc];
	do
		dslc = slc(i+1)-slc(i);
		V1 = shock_design(rlc(i,:),slc(i),clc,s0lc,r_uc).*dslc;
		V2 = shock_design(rlc(i,:)+V1./2,slc(i)+dslc/2,clc,s0lc,r_uc).*dslc;
		V3 = shock_design(rlc(i,:)+V2./2,slc(i)+dslc/2,clc,s0lc,r_uc).*dslc;
		V4 = shock_design(rlc(i,:)+V3,slc(i)+dslc,clc,s0lc,r_uc).*dslc;
		rlc(i+1,:) = (V1+2.*V2+2.*V3+V4)./6+rlc(i,:);
		i++;
	until (slc(i) == 0 || i>=10000)
	res =rlc(rows(rlc),2)/rlc(rows(rlc),1)
end
	
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

function res = cone_solve(V0,cone_theta,kappa,beta)
	theta = linspace(beta,cone_theta,100)';
	Vd = lsode(@(V,theta)taylor_maccol(V,theta,kappa),V0,theta);
	rad2deg(atan(Vd(rows(theta),2)/Vd(rows(theta),1)));
	(180-rad2deg(cone_theta));
	res = rad2deg(atan(Vd(rows(theta),2)/Vd(rows(theta),1)))-(180-rad2deg(cone_theta));
end

function [z_out,y_out]=streamline(r0,theta,Vr,Vs,zte)
	r(1) = r0;
	for i= 1:rows(Vr)-1
		dtheta(i) = -(theta(i+1)-theta(i));
		a(i) = -atan(Vr(i)/Vs(i))+pi/2;
		dr(i) = (sin(a(i))/sin(dtheta(i)+a(i))-1)*r(i);
		r(i+1) = r(i)+dr(i);
	end

	z= r'.*cos(theta);
	y=r'.*sin(theta);
	
	%長さをそろえる
	z_out = z;
	y_out = y;
	
end
	
%勾配法により中心点を求める
function res = circle_solve_x(cx,x_el,y_el,psi,R)
	res = R^2-((x_el-cx)^2+(y_el-(tan(psi)*(cx-x_el)+y_el))^2);
end
function res = circle_solve_y(cy,x_el,y_el,psi,R)
	res = R^2-((x_el-((cy-y_el)/tan(psi)+x_el))^2+(y_el-cy)^2);
end