1;
clear all


function y = machplot(H,am,x0,y0,x)
	y = -tan(am-H+deg2rad(0)).*(x-x0)+y0;
end




%prandtl meyerによるインテーク設計
M_mainflow =input("Mach mainflow: ");
%流入角
aoa = -input('aoa(deg) ' )
disp('----cf----')
beta_cf=fsolve(@(beta)tbm_function(deg2rad(-aoa),M_mainflow,beta,1.4),0.0001);
rad2deg(beta_cf)
[Min_cf p02_p01_cf p2_p1_cf]=pmvalue(M_mainflow,beta_cf,1.4)
disp('----------')
%機首予圧縮
H0 = -input("H 0(deg): ");

Min =input("Mach in: ");

if aoa-H0 != 0
	beta_nose = fsolve(@(beta)tbm_function(deg2rad(-H0+aoa),Min,beta,1.4),0.0001);
	[Min_nose p02_p01_nose p2_p1_nose]=pmvalue(Min,beta_nose,1.4);
	beta_nose= beta_nose-deg2rad(aoa)
else;
	beta_nose = atan(1./sqrt(Min.^2-1))-deg2rad(aoa);
	Min_nose = Min;
	p02_p01_nose = 1; 
	p2_p1_nose  = 1;
end
Min_nose
Mout =input("Mach out: ");
Tmc = input("input Tmc ");
%Mout_list = linspace(1.01,5.9,100)';




hin =3660 %73; %input("inlet hight:")
p0 = 5531;
n_div =10000;
%k=0.985

%~ Min =6 %input("Mach in: ");
%~ Mout =3.95%input("Mach out: ");
%~ H0 = -input("H 0(deg): ")

%~ hin = 3660; %input("inlet hight:")
%~ ht =1140 %input("throat hight:")
%~ p0 = 5531
%~ n_div =30;
% k = 0.96

%k = fsolve(@(k)pminlet_k_solve(Min,Mout,n_div,H0,hin,p0,k),k)


theta0 = prandtl_meyer_function(Min_nose,1.4)-prandtl_meyer_function(Mout,1.4); %deg
%dtheta = -theta0/(n_div)*k
%phi0 = 170
%res = pminlet_solve(Min,Mout,n_div,H0,phi0)
%phi0 = fsolve(@(phi0)pminlet_solve(Min,Mout,n_div,H0,phi0,deg2rad(dtheta)),179)


M = linspace(Min_nose,Mout,n_div)';
am =  atan(1./sqrt(M.^2-1));
phi0 =  180 - rad2deg(am(1,1)) + H0;
theta(1,1) =deg2rad(prandtl_meyer_function(M(2,1),1.4)-prandtl_meyer_function(M(1,1),1.4));
r(1,1) = 1;
p(1,1) = p0*p2_p1_nose*p2_p1_cf;
phi(1,1) = deg2rad(phi0);
for i = 2:n_div
	theta(i,1) = deg2rad(prandtl_meyer_function(M(i,1),1.4) - prandtl_meyer_function(M(i-1,1),1.4));
	dp(i,1) = -theta(i,1)*1.4*M(i,1)^2/sqrt(M(i,1)^2-1)*p(i-1,1);
	p(i,1) = p(i-1,1)+dp(i,1);
	dr(i,1) = (sin(am(i-1,1)+theta(i-1,1))/sin(am(i,1))-1)*r(i-1,1);
	r(i,1) = r(i-1,1)+dr(i,1);
	dphi(i,1) = am(i,1)-(am(i-1,1)+theta(i-1,1));
	phi(i,1) = phi(i-1,1)-dphi(i,1);
end
x = r.*cos(phi);
y = r.*sin(phi);

H(1,1) = deg2rad(H0);
for i = 2:n_div
	H(i,1) = atan((y(i)-y(i-1))/(x(i)-x(i-1)));
end

%インテークの相似拡大
x_intake = hin./max(y) .*x;
y_intake = hin./max(y) .*y;
r0 = sqrt(x_intake(1,1)^2+y_intake(1,1)^2);

%接続部は5次関数

%衝撃波角計算
tan_slast = tan(H(rows(H)));%(y_intake(rows(y),1)- y_intake(rows(y)-1,1))/(x_intake(rows(x),1)- x_intake(rows(x)-1,1));
beta = fsolve(@(beta)tbm_function(-atan(tan_slast),M(rows(M),1),beta,1.4),0.001);
shock_angle = beta+H(rows(H),1);
[mach2 p02_p01 p2_p1] = pmvalue(M(rows(M),1),beta,1.4);
MCR = 1
TPR = p02_p01*p02_p01_nose*p02_p01_cf
mach2 
p_pin = p2_p1*p(rows(p),1)/p0
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
At_Ain = ys/Ain
beta =rad2deg(beta)
theta0
shock_theta =- rad2deg(atan(tan_slast))

%ノーズ位置計算
x_nose = (-x_intake(1,1)*tan(H(1,1))+y_intake(1,1))/(-tan(H(1,1))-tan(beta_nose));
y_nose = tan(-beta_nose)*x_nose;


%~ %接続部計算
%~ A = [xs^3     , xs^2 , xs, 1;
       %~ xo^3     , xo^2, xo,  1;
       %~ 3*xs^2 ,2*xs  ,1  ,  0;
       %~ 3*xo^2 ,2*xo  ,1  ,  0];
%~ B = [ys;
       %~ yo;
       %~ 0;
      %~ (y_intake(rows(y),1)- y_intake(rows(y)-1,1))/(x_intake(rows(x),1)- x_intake(rows(x)-1,1))];
      
%~ par = A\B;
%~ x_sp = linspace(xo,xs,500)';
%~ y_sp = par(1) .*x_sp .^ 3 + par(2) .* x_sp.^ 2 + par(3) * x_sp.^ 1 + par(4);



%マッハ波の描画


plot(x_intake,y_intake,'r')
hold on
plot([x_nose,x_intake(1,1)],[y_nose,y_intake(1,1)],'k')
plot([x_nose,x_nose-hin],[y_nose,tan(deg2rad(aoa))*(-hin)+y_nose],'c')
plot([0,x_nose],[0,y_nose],'m')
plot([0,x_sholder],[0,ht],"r")
plot(x_sp,y_sp)
plot([xs xs+ht],[ys,ys])
plot([0 xs+ht],[0,0])

printf("shock angle error: %f \n",shock_theta-theta0+H0)


for i = 2:1000:n_div-1
	x_plot= linspace(x_intake(i,1),0,100)';
	y_plot = machplot(H(i,1),am(i,1),x_intake(i,1),y_intake(i,1),x_plot);
	plot(x_plot,y_plot,'g')
end
	i = n_div;
	x_plot= linspace(x_intake(i,1),0,100)';
	y_plot = machplot(H(i,1),am(i,1),x_intake(i,1),y_intake(i,1),x_plot);
	plot(x_plot,y_plot,'g')

axis "equal"
ylim([0,max(y_intake)*1.5])
hold off
	
%出力ファイルの用意
x_out = linspace(min(x_intake),max(x_intake),200)';
y_out =spline([x_intake],[y_intake],x_out);
dlmwrite("pminlet.txt",[[0;xs;x_nose;x_out],[0;ys;y_nose;y_out]],delimiter = " ")


%推力計算

[F Isp] = scram_thrust(MCR,TPR,p_pin,Tmc,M_mainflow,216,5531,0.088,max(y_nose)/1000*2.02,ys/1000*2.02,3.66*2.02)
		

%------他設計点解析



		