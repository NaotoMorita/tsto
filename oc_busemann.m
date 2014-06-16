1;
clear all
function dV =  taylor_maccoll(V,theta,kappa)
	dV(1) = V(2);
	x=(kappa-1)/2* (1-V(1)^2-V(2)^2);
	dV(2) = (V(2)^2*V(1)-2*x*V(1)-x*V(2)/tan(theta)) / (x-V(2)^2);
end

function res = cone_solve(V0,cone_theta,kappa,beta)
	theta = linspace(beta,cone_theta,100)';
	Vd = lsode(@(V,theta)taylor_maccol(V,theta,kappa),V0,theta);
	rad2deg(atan(Vd(rows(theta),2)/Vd(rows(theta),1)))
	(180-rad2deg(cone_theta))
	res = rad2deg(atan(Vd(rows(theta),2)/Vd(rows(theta),1)))-(180-rad2deg(cone_theta))
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
kappa = 1.4

beta = deg2rad(34.3797988108415)
M1 = 2.5


delta = atan(2*cot(beta)*((M1^2*sin(beta)^2-1)/(M1^2*(kappa+cos(2*beta))+2)));
Mn1 = M1*sin(beta);
Mn2 = sqrt((Mn1^2+(2/(kappa-1)))/(2*kappa*Mn1^2/(kappa-1)-1));
M2 = Mn2/sin(beta-delta);

Vd = ( 2/(kappa-1)/M1^2+1)^(-1/2);  
Vdr0 = Vd*cos(beta);
Vdt0 = -Vd*sin(beta);
Vrs0=[Vdr0,Vdt0];

s(1) = beta-delta;
ds = deg2rad(0.05);
Vrs(1,1:2) = Vrs0;
Vxy(1,1) = Vrs0(1,1) *cos(s(1))- Vrs0(1,2) *sin(s(1));
Vxy(1,2) = Vrs0(1,1) *sin(s(1))+ Vrs0(1,2) *cos(s(1));
r(1)=1;
z(1) = r(1)*cos(s(1));
x(1) = r(1)*sin(s(1));
Vabs(1) = sqrt(Vrs(1,1)^2+Vrs(1,2)^2)
mach(1) = M1

%4次rungekuttaによるTM方程式の数値解
i=1;
do
	s(i+1) =  s(i)+ds;
	V1 = taylor_maccoll(Vrs(i,:),s(i),1.4).*ds;
	V2 = taylor_maccoll(Vrs(i,:)+V1./2,s(i)+ds/2,1.4).*ds;
	V3 = taylor_maccoll(Vrs(i,:)+V2./2,s(i)+ds/2,1.4).*ds;
	V4 = taylor_maccoll(Vrs(i,:)+V3,s(i)+ds,1.4).*ds;
	Vrs(i+1,:) = (V1+2.*V2+2.*V3+V4)./6+Vrs(i,:);
	
	Vabs(i+1) = sqrt(Vrs(i+1,1)^2+Vrs(i+1,2)^2);
	mach(i+1) = sqrt(2/((1/Vabs(i+1)^2-1)*(kappa-1)));
	Vxy(i+1,1) = Vrs(i+1,1) *cos(s(i+1)) - Vrs(i+1,2) *sin(s(i+1));
	Vxy(i+1,2) = Vrs(i+1,1) *sin(s(i+1)) + Vrs(i+1,2) *cos(s(i+1));
	
	%streamlineも作る
	if -atan(Vrs(i+1,2)/Vrs(i+1,1)) >0
		a = -atan(Vrs(i+1,2)/Vrs(i+1,1));
	else
		a = -atan(Vrs(i+1,2)/Vrs(i+1,1))+pi;
	end
	dr = (sin(a)/sin(ds+a)-1)*r(i);
	r(i+1) = r(i)+dr;
	z(i+1) = r(i+1)*cos(s(i+1));
	x(i+1) = r(i+1)*sin(s(i+1));
i++
until (Vxy(i,2) >= -0.001 || i>=10000)
i--;
z -= z(1,1);
x -= x(1,1);

%楕円スロートの曲率半径を計算する
%楕円

a_el = 0.914;
b_el = 1;



%入口矩形インテーク
x_rec = 0.5;
y_rec = 0.75%b_el*sqrt(1-(x_rec/a_el)^2);

%境目を見つける
function res = border_solve(x_bo,a,b,x_rec,y_rec)
	y_bo = b*sqrt(1-(x_bo/a)^2);
	phi = atan(y_bo/x_bo);
	k = a*b/(a^2*sin(phi)^2+b^2*cos(phi)^2)^(3/2);
	R = 1/k;
	psi = atan(a*y_bo/(b*x_bo));
	cx = fsolve(@(cx)circle_solve_x(cx,x_bo,y_bo,psi,R),-10);
	cy = tan(phi)*(cx-x_bo)+y_bo;
	if cx == x_bo
		cy = fsolve(@(cy)circle_solve_y(cy,x_bo,y_bo,psi,R),-10);
		cx = (cy-y_bo)/tan(psi)+x_bo;
	end
	xstt = x_rec;
	ystt =  tan(psi)*(x_rec-x_bo)+y_bo;
	res = ystt-y_rec
end

x_bo = fsolve(@(x_bo)border_solve(x_bo,a_el,b_el,x_rec,y_rec),0.5)
y_bo = b_el*sqrt(1-(x_bo/a_el)^2);
phi_bo = atan(sqrt(a_el^2-x_bo^2)/x_bo)
phi = linspace(0,phi_bo,20)';
phi(rows(phi),:)=[];
phi(20:34,1) = linspace(phi_bo,pi/2,15)';

c_el = sqrt(abs(a_el.^2-b_el.^2))
r_el = b_el.^2./(a_el+c_el.*cos(phi));

x_el =a_el.*cos(phi);
y_el =b_el.*sin(phi);
[buff i]= sort(x_el, 'descend');



R_real = sqrt(x_el.^2+y_el.^2);
z_el(1:rows(y_el),1) = 0;
%曲率
k_el = a_el.*b_el./(a_el.^2.*sin(phi).^2+b_el.^2.*cos(phi).^2).^(3/2);
R_el = 1./k_el
printf('mach: %f',mach(length(mach)))


for j = 1:rows(x_el)
	%局所円錐衝撃波面の傾きψ
	psi(j,1) = atan(a_el*y_el(j)/(b_el*x_el(j)));
	%中心点
	cx(j) = fsolve(@(cx)circle_solve_x(cx,x_el(j),y_el(j),psi(j),R_el(j)),-10);
	cy(j) = tan(psi(j))*(cx(j)-x_el(j))+y_el(j);
	if cx(j) == x_el(j)
		cy(j) = fsolve(@(cy)circle_solve_y(cy,x_el(j),y_el(j),psi(j),R_el(j)),-10);
		cx(j) = (cy(j)-y_el(j))/tan(psi(j))+x_el(j);
	end
	cz(i) = R_el/sin(beta-delta);

end


%OCBIを形作る
	for j = 1:rows(x_el)
		if x_el(j)>x_bo
			xstt(j) = x_rec;
			ystt(j) = tan(psi(j))*(x_rec-x_el(j))+y_el(j);
			r_rec(j) = sqrt((xstt(j)-cx(j))^2+(ystt(j)-cy(j))^2);
		else
			%ystt(j) = y_el(j);
			%xstt(j) = x_el(j);
			ystt(j) = y_rec;
			xstt(j) = (y_rec-y_el(j))/tan(psi(j))+x_el(j);
			r_rec(j) = sqrt((xstt(j)-cx(j))^2+(ystt(j)-cy(j))^2);
		end
		zstt(j) = (R_el(j)-r_rec(j))/tan(beta-delta);
		rcone(j) =r_rec(j)/sin(beta-delta);
		x_oc(:,j) = xstt(j)+rcone(j).*(x'.*cos(psi(j))); 
		y_oc(:,j) = ystt(j)+rcone(j).*(x'.*sin(psi(j))); 
		z_oc(:,j) = -zstt(j)+rcone(j).*R_el(j).*(z');
		%x_oc(i,j) = x_el(j)+R_el(j)*(x(i)*cos(phi(j))); 
		%y_oc(i,j) = y_el(j)+R_el(j)*(x(i)*sin(phi(j))); 
		%z_oc(i,j) = z_el(j)+R_el(j)*(z(i));
	end


figure(1)
nowwork=pwd;
cd OCBI
for j = 1:rows(x_el)
	plot3(x_oc(:,j),y_oc(:,j),z_oc(:,j))
	hold on
	plot3(cx(j),cy(j),0,'go')
	plot3([xstt(j),cx(j)],[ystt(j),cy(j)],[0,0],'-go')
	
	z_out = linspace(z_oc(1,j),z_oc(rows(x_oc),j),200)';
	x_out = spline(z_oc(:,j),x_oc(:,j),z_out);
	y_out = spline(z_oc(:,j),y_oc(:,j),z_out);
	if min(x_oc(:,j)) <= 0.0001
		1
		dlmwrite(strcat('stream',mat2str(j),'.sldcrv'),[zeros(rows(x_out),1),y_out,z_out].*1000);
	else
		dlmwrite(strcat('stream',mat2str(j),'.sldcrv'),[x_out,y_out,z_out].*1000);
	end
	plot3(-x_oc(:,j),y_oc(:,j),z_oc(:,j))
	
	
end
cd(nowwork)
hold off
axis 'equal'
figure(2)
	plot(xstt,ystt)
	hold off
	axis 'equal'
	
%衝撃波形状
figure(4)
for j = 1:rows(x_el)
	plot3([cx(j),xstt(j)],[cy(j),ystt(j)],[cz(j),zstt(j)],'r')
	plot3(-[cx(j),xstt(j)],[cy(j),ystt(j)],[cz(j),zstt(j)],'r')
	hold on
	plot3([x_el(j),xstt(j)],[y_el(j),ystt(j)],[0,zstt(j)])
	plot3(-[x_el(j),xstt(j)],[y_el(j),ystt(j)],[0,zstt(j)])
	
end
axis 'equal'
view(210,45)
grid on
hold off