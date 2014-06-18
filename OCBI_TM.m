
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