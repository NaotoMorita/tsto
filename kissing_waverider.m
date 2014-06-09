1;
clear all
function dV =  taylor_maccol(V,theta,kappa)
	dV(1) = V(2);
	x=(kappa-1)/2* (1-V(1)^2-V(2)^2);
	dV(2) = (V(2)^2*V(1)-2*x*V(1)-x*V(2)/tan(theta)) / (x-V(2)^2);
end

function res = cone_solve(V0,cone_theta,kappa,beta)
	theta = linspace(beta,cone_theta,100)';
	Vd = lsode(@(V,theta)taylor_maccol(V,theta,kappa),V0,theta);
	res = Vd(rows(theta),2);
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
	y=-r'.*sin(theta);
	y -= y(1,1);
	z -= z(1,1);
	
	%長さをそろえる
	z_out = linspace(0,zte,500);
	y_out = interp1(z,y,z_out);
	
end

%三点から円の中心と半径を求める
function [ac,r] = radius_center(r1,r2,r3)
	A = [r1(1),r1(2),1;
	       r2(1),r2(2),1;
	       r3(1),r3(2),1];
	B=[-r1(1)^2-r1(2)^2;
	     -r2(1)^2-r2(2)^2;
	     -r3(1)^2-r3(2)^2];
	lmn = A\B;
	l = lmn(1);
	m= lmn(2);
	n = lmn(3);
	
	ac(1) = -l/2;
	ac(2) = -m/2;
	r = sqrt((l/2)^2+(m/2)^2-n);;
end


kappa = 1.4
beta = deg2rad(13)
M1 = 6

delta = atan(2*cot(beta)*((M1^2*sin(beta)^2-1)/(M1^2*(kappa+cos(2*beta))+2)));
Mn1 = M1*sin(beta);
Mn2 = sqrt((Mn1^2+(2/(kappa-1)))/(2*kappa*Mn1^2/(kappa-1)-1));
M2 = Mn2/sin(beta-delta);

Vd = ( 2/(kappa-1)/M2^2+1)^(-1/2);  
Vdr0 = Vd*cos(beta-delta);
Vdt0 =  -Vd*sin(beta-delta);
V0=[Vdr0,Vdt0];
cone_theta = fsolve(@(cone_theta)cone_solve(V0,cone_theta,kappa,beta),beta-0.01);
theta = linspace(beta,cone_theta,100)';
[Vd] = lsode(@(V,theta)taylor_maccol(V,theta,kappa),V0,theta);
Vr = Vd(:,1);
Vs = Vd(:,2);

Vx = Vr.*cos(theta)-Vs.*sin(theta);
Vy = Vr.*sin(theta)+Vs.*cos(theta);

%shockline
Xshock = 0.13
Yshock = 0.25
Xbody = 0.30
Ybody = 0.65
Xinlet = 0.17
Yinlet = 0.65
shock_dd = -1

n_line = 20;
n_curve = 502;
x_sl = linspace(0,Xshock,n_line)';
y_sl(1:n_line,1) = 0;

%二次関数を作成
A=[Xshock^3,Xshock^2,Xshock,1;
     1,          1            ,1        ,1;
     3*Xshock^2,2*Xshock,1        ,0;
     6*Xshock,   2,             0,        0]
B=[0;
     Yshock;
     0;
     0]
pol = A\B;
th = linspace(0,pi,n_curve)';
x_sc(:,1) = (1-(cos(th)./2+0.5).^2).*(1-Xshock)+Xshock;
y_sc(:,1) = pol(1).*x_sc(:,1).^3+pol(2).*x_sc(:,1).^2+pol(3).*x_sc(:,1)+pol(4);
%曲率から中心を求める
for i= 1:n_curve-1
	if i==1
		[ac(i,:),r(i)]=radius_center([x_sc(1,1),y_sc(1,1)],[x_sc(2,1),y_sc(2,1)],[x_sc(3,1),y_sc(3,1)]);;
	else 
		[ac(i,:),r(i)]=radius_center([x_sc(i-1,1),y_sc(i-1,1)],[x_sc(i,1),y_sc(i,1)],[x_sc(i+1,1),y_sc(i+1,1)]);;
	end
	%面の傾きを求めておく
	phi(i) = atan((y_sc(i,1)-ac(i,2))/(ac(i,1)-x_sc(i,1)));
	
	%前縁と曲率半径の接する線の判別%とりあえず前縁は直線とする
	Yinlet_line(i) = tan(-phi(i))*(Xinlet-ac(i,1))+ac(i,2);
	if Yinlet_line(i) <= Yinlet
		gradbuff = (Yinlet-Yinlet)/(Xinlet-0);
		X_LE(i) = (tan(-phi(i))*ac(i,1)-ac(i,2)+Yinlet)/(tan(-phi(i))-gradbuff);
		Y_LE(i) = Yinlet;
		no(i) = 1; 
	else
		Ybody_line(i) = tan(-phi(i))*(Xbody-ac(i,1))+ac(i,2);
		if Ybody_line(i)<=Ybody
			gradbuff = (Ybody-Yinlet)/(Xbody-Xinlet);
			X_LE(i) = (tan(-phi(i))*ac(i,1)-ac(i,2)-gradbuff*Xbody+Ybody)/(tan(-phi(i))-gradbuff);
			Y_LE(i) = tan(-phi(i))*(X_LE(i)-ac(i,1))+ac(i,2);
			no(i) = 2; 
		else
			gradbuff = (Yshock-Ybody)/(1-Xbody);
			X_LE(i) = (tan(-phi(i))*ac(i,1)-ac(i,2)-gradbuff*1+Yshock)/(tan(-phi(i))-gradbuff);
			Y_LE(i) = tan(-phi(i))*(X_LE(i)-ac(i,1))+ac(i,2);
			no(i) = 3; 
		end
	end
	
	
	%ZLE前縁の奥行
	Z_LE(i) = cot(beta)*(sqrt((X_LE(i)-x_sc(i))^2+(Y_LE(i)-y_sc(i))^2));	
	
	%とりあえず、後縁形状を調べる
	Z_cone(i) = cot(beta)*(sqrt((ac(i,1)-x_sc(i))^2+(ac(i,2)-y_sc(i))^2));
	R_stream(i) = (sqrt((ac(i,1)-X_LE(i))^2+(ac(i,2)-Y_LE(i))^2))/sin(beta);
	[Z_str,Y_str] = (streamline(R_stream(i),theta,Vr,Vs,Z_LE(i)));
	Y_TE(i) = min(Y_str.*sin(phi(i)))+Y_LE(i);
	Y_stream(:,i) = Y_str.*sin(phi(i))+Y_LE(i);
	X_stream(:,i) = -Y_str.*cos(phi(i)) +X_LE(i);
	Z_stream(:,i) = Z_LE(i)-Z_str;
	X_TE(i) = (Y_TE(i)-ac(i,2))/tan(-phi(i))+ac(i,1);
end
Y_TE(1) = Y_TE(2);
X_TE(1) = X_TE(2);
Z_LE(1) = Z_LE(2);
Y_LE(1) = Y_LE(2);
X_LE(1) = X_LE(2);
Y_stream(:,1) = Y_stream(:,2);
X_stream(:,1) = X_stream(:,2);
Z_stream(:,1) = Z_stream(:,2);
figure(1)
plot([0,Xinlet],[Yinlet,Yinlet])
hold on
plot([Xinlet,Xbody],[Yinlet,Ybody],'r--')
plot([Xbody,1],[Ybody,Yshock],'r--')
plot(X_LE,Y_LE)
plot([x_sl;x_sc],[y_sl;y_sc],'r')
%plot(ac(:,1),ac(:,2))
plot(X_TE,Y_TE,'-o')
plot([0,X_TE(1,1)],[Y_TE(1,1),Y_TE(1,1)],'-o')
axis 'equal'
hold off

figure(2)
%~ hold off
%~ for i = 2:10:n_curve-1
	%~ plot3(X_stream(:,i),Y_stream(:,i),Z_stream(:,i),'r')
	%~ plot3(-X_stream(:,i),Y_stream(:,i),Z_stream(:,i),'r')
	%~ plot3([X_LE(i),X_LE(i)],[Y_LE(i),Y_LE(i)],[0,Z_LE(i)])
	%~ plot3([-X_LE(i),-X_LE(i)],[Y_LE(i),Y_LE(i)],[0,Z_LE(i)])
	%~ filename = strcat("waveride_sf",mat2str(i),'.sldcrv')
	%~ nowwork = pwd;
	%~ cd("waverider")
	%~ upsf(1:150,3) = linspace(0,Z_LE(i),150)';
	%~ upsf(1:150,2)=Y_LE(i);
	%~ upsf(1:150,1)=X_LE(i);
	%~ upsf(150,:)=[];
	%~ endsf (:,1) = linspace(X_TE(i),X_LE(i),150)';
	%~ endsf (:,2) = linspace(Y_TE(i),Y_LE(i),150)';
	%~ endsf(:,3)=0;
	%~ endsf(150,:)=[];
	%~ losf = [X_stream(:,i),Y_stream(:,i),Z_stream(:,i)];
	%~ dlmwrite(filename,[losf])
	%~ cd(nowwork)
	%~ clear upsf endsf losf
	%~ hold on 
	%~ i
%~ end
%~ cd("waverider")
%~ i=n_curve-1
%~ filename = strcat("waveride_sf",mat2str(i),'.sldcrv')
%~ losf = [X_stream(:,i),Y_stream(:,i),Z_stream(:,i)];
%~ dlmwrite(filename,[losf])
%~ cd(nowwork)

%~ %直線部分を描画
%~ for i = 1:1:n_line
	%~ [Z_str,Y_str] = streamline(R_stream(2),theta,Vr,Vs,Z_LE(2));
	
	%~ Y_stl(:,i) = Y_str+Yinlet;
	%~ X_stl(1:rows(Y_stl(:,i)),i) = x_sl(i);
	%~ Z_stl(1:rows(Y_stl(:,i)),i) = Z_LE(2)-Z_str;
	%~ plot3(X_stl(:,i),Y_stl(:,i),Z_stl(:,i),'r')
	%~ plot3(-X_stl(:,i),Y_stl(:,i),Z_stl(:,i),'r')
	%~ plot3([x_sl(i),x_sl(i)],[Yinlet,Yinlet],[0,Z_LE(2)])
	%~ plot3([-x_sl(i),-x_sl(i)],[Yinlet,Yinlet],[0,Z_LE(2)])
%~ end
%~ axis"equal"
%~ hold off

%断面作成
n = 1
j = 1
k = 1
wr_splot(1) = 0;
do
	n
	if Z_LE(n) == Z_LE(n+1)
		n++;
	else
		if wr_splot(1) == 0;
			wr_splot(1) = j+1
		end
		%~ %直線部分の補間
		%~ slice_wave{j}(1:n_line,1) = linspace(0,interp1(Z_stream(:,1),X_stream(:,1),Z_LE(n),'extrap'),n_line)';
		%~ slice_wave{j}(1:n_line,2) = interp1(Z_stream(:,1),Y_stream(:,1),Z_LE(n),'extrap');
		%~ slice_wave{j}(1:n_line,3) = Z_LE(n);
		
			
		for i = 1:n
			slice_wave{j}(i,1) = interp1(Z_stream(:,i),X_stream(:,i),Z_LE(n),'extrap');
			slice_wave{j}(i,2) = interp1(Z_stream(:,i),Y_stream(:,i),Z_LE(n),'extrap');
			slice_wave{j}(i,3) = Z_LE(n);
		end
		slice_wave{j}(1,:)=[];
		%~ clear slice_wavebuff
		for i = 2:n
			slice_LE{j}(i,1) =X_LE(i);
			slice_LE{j}(i,2) =Y_LE(i);
			slice_LE{j}(i,3) =Z_LE(n);
		end
		%~ slice_wavebuff(1,:) = [];
		
		%~ slice_wave{j}(n_line+n-2:n_line+n+rows(slice_wavebuff)-3,:)=flipud(slice_wavebuff);
		%~ slice_wave{j}(rows(slice_wave{i}),:)=[];
		%~ %直線部分の補間
		%~ slice_wave{j}(n_line+n+rows(slice_wavebuff)-3:n_line+n+rows(slice_wavebuff)+n_line-4,1) = linspace(X_LE(1),0,n_line)';
		%~ slice_wave{j}(n_line+n+rows(slice_wavebuff)-3:n_line+n+rows(slice_wavebuff)+n_line-4,2) = Yinlet;
		%~ slice_wave{j}(n_line+n+rows(slice_wavebuff)-3:n_line+n+rows(slice_wavebuff)+n_line-4,3) = Z_LE(n);
		%~ slice_wave{j}(rows(slice_wave{j})+1:rows(slice_wave{j})*2,:) =flipud([-slice_wave{j}(:,1),slice_wave{j}(:,2),slice_wave{j}(:,3)]);
		%~ slice_wave{j}(rows(slice_wave{j})/2,:) =[];
		
		%間引きを行う
		dz = Z_LE(1) / 4 * (n/rows(Z_LE'))^1.7;
		Z_LE(wr_splot(k))-Z_LE(n)
		if  Z_LE(wr_splot(k))> Z_LE(n)+dz || n == rows(Z_LE') || n == rows(Z_LE')-1 
			k++
			wr_splot(k) = j;
			
		end
				n++;
				j++;
	end
until (n== rows(Z_LE'))
j--;
n--;
nowwork = pwd;
cd("waverider")	
for n = 1:length(wr_splot)

	plot3(slice_wave{wr_splot(n)}(:,1),slice_wave{wr_splot(n)}(:,2),slice_wave{wr_splot(n)}(:,3))
	if n == 1
		i=-1
		do
			i++
		until (max(slice_wave{wr_splot(n)+i}(:,2)) - min(slice_wave{wr_splot(n)+i}(:,2))>10^(-3) || wr_splot(n)+i == wr_splot(n+1))
		filename = strcat("waveride_sf",mat2str(n),'.sldcrv')
		dlmwrite(filename,slice_wave{wr_splot(n)+i}.*1000)
	elseif n== length(wr_splot)
		filename = strcat("waveride_sf",mat2str(n),'.sldcrv')
		dlmwrite(filename,[slice_wave{wr_splot(n)}(:,1:2),zeros(rows(slice_wave{wr_splot(n)}(:,1:2)),1)].*1000)
	else
		filename = strcat("waveride_sf",mat2str(n),'.sldcrv')
		dlmwrite(filename,slice_wave{wr_splot(n)}.*1000)
	end
	hold on
end

cd(nowwork)
axis 'equal'
hold off


