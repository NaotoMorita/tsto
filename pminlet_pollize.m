%PMinletを多角形化し、性能比較するコード
%多角形化数
n_pol = 2

clear x_pol y_pol M_pol beta_pol p02_p01_pol p2p1_pol pan_angle

%等エントロピー圧縮部分の多角形化
for i = 2:n_pol+2
	x_pol(i,1)= -(x_intake(1,1)-x_intake(rows(x_intake),1))/(n_pol+1)*(i-1)+x_intake(1,1);
end
y_pol(1:rows(x_pol)-1,1) = interp1(x_intake(:,1),y_intake(:,1),x_pol(1:rows(x_pol)-1,1));
x_pol(1,1) = x_intake(1,1);
x_pol(rows(x_pol),1)=x_intake(rows(x_intake),1);
y_pol(1,1) = y_intake(1,1);
y_pol(rows(x_pol),1)=y_intake(rows(y_intake),1);

%転向角と衝撃波計算
%パネル角度計算
for i = 2:n_pol+2
	pan_angle(i-1) = atan((y_pol(i,1)-y_pol(i-1,1))/(x_pol(i,1)-x_pol(i-1,1)));
end	

M_pol(1) = Min_nose
for i = 1:n_pol+2
	if i == 1
		beta_pol(i)=fsolve(@(beta)tbm_function(+deg2rad(H0)-pan_angle(i),M_pol(i),beta,1.4),0.0001);
	elseif i== n_pol+2
		beta_pol(i)=fsolve(@(beta)tbm_function(-H(rows(H))+pan_angle(length(pan_angle)),M_pol(i),beta,1.4),0.0001);
	else
		beta_pol(i)=fsolve(@(beta)tbm_function(-pan_angle(i)+pan_angle(i-1),M_pol(i),beta,1.4),0.0001);
	end
	[M_pol(i+1) p02_p01_pol(i) p2_p1_pol(i)]=pmvalue(M_pol(i),beta_pol(i),1.4);
end

%終端衝撃波
betaT_pol = fsolve(@(beta)tbm_function(-atan(tan_slast),M_pol(length(M_pol)),beta,1.4),0.001);
[mach2_pol p02_p01T_pol p2_p1T_pol] = pmvalue(M_pol(length(M_pol)),betaT_pol,1.4);

TPR_pol = p02_p01_cf*p02_p01_nose*p02_p01T_pol;
p_pin_pol = p2_p1T_pol*p2_p1_nose*p2_p1_cf;
for i = 1:rows(n_pol+2)
	TPR_pol *= p02_p01_pol(i);
	p_pin_pol *= p2_p1_pol(i);
end
mach2_pol

%以下描画
plot(x_intake,y_intake,'r')
hold on
plot([x_nose,x_intake(1,1)],[y_nose,y_intake(1,1)],'k')
plot([x_nose,x_nose-hin],[y_nose,tan(deg2rad(aoa))*(-hin)+y_nose],'c')
plot([0,x_nose],[0,y_nose],'m')
plot([0,x_sholder],[0,ht],"r")
plot(x_sp,y_sp)
plot([xs xs+ht],[ys,ys])
plot([0 xs+ht],[0,0])

for i = 2:n_pol+2
	x_plot= linspace(x_pol(i,1),0,100)';
	y_plot = machplot(pan_angle(i-1),beta_pol(i),x_pol(i,1),y_pol(i,1),x_plot);
	plot(x_plot,y_plot,'g')
end
	i = 1;
	x_plot= linspace(x_pol(i,1),0,100)';
	y_plot = machplot(deg2rad(H0),beta_pol(i),x_pol(i,1),y_pol(i,1),x_plot);
	plot(x_plot,y_plot,'g')
	

axis 'equal'
hold on
plot(x_pol,y_pol)
hold off

[F Isp] = scram_thrust(MCR,TPR_pol,p_pin_pol,Tmc,M_mainflow,216,5531,0.088,max(y_nose)/1000*2.02,ys/1000*2.02,3.66*2.02)