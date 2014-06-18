
%スロート&衝撃波設計
%入口矩形インテーク
x_rec = 0.5;
y_rec = 0.75;

%上部円、下曲率の衝撃波形状
y_uc = fsolve(@(a)shock_upper(a,x_rec,y_rec),0.01); 
r_uc = sqrt(x_rec^2+(y_rec-y_uc)^2);
x_cone = linspace(0,x_rec,20)';
y_cone = y_uc+(sqrt(r_uc.^2-x_cone.^2));
x_cone(rows(x_cone),:)=[];
y_cone(rows(y_cone),:)=[];
%4次rungekuttaによるTM方程式の数値解
i=1;
r(1) = sqrt(y_rec^2+x_rec^2);
s0lc = atan(y_rec/x_rec);
slc= linspace(atan(y_rec/x_rec),0,200)';
rlc(1,:) = [r(1),(y_uc*r(1)*cos(s0lc))/(r(1)-y_uc*sin(s0lc)),r_uc]
clc = fsolve(@(clc)clc_solve(clc,x_rec,y_rec,r_uc,y_uc),-0.5)
do
	dslc = slc(i+1) -  slc(i);
	V1 = shock_design(rlc(i,:),slc(i),clc,s0lc,r_uc).*dslc;
	V2 = shock_design(rlc(i,:)+V1./2,slc(i)+dslc/2,clc,s0lc,r_uc).*dslc;
	V3 = shock_design(rlc(i,:)+V2./2,slc(i)+dslc/2,clc,s0lc,r_uc).*dslc;
	V4 = shock_design(rlc(i,:)+V3,slc(i)+dslc,clc,s0lc,r_uc).*dslc;
	rlc(i+1,:) = (V1+2.*V2+2.*V3+V4)./6+rlc(i,:);
	i++;
until (slc(i) <= 0 || i>=10000)
res_tre =rlc(rows(rlc),2)/rlc(rows(rlc),1)
slc(i)=[];
rlc(i,:) =[];

slcout = linspace(atan(y_rec/x_rec),0,20)';
rlcout = interp1(slc',rlc,slcout,'extrap');

xlc = rlcout(:,1).*cos(slcout)
ylc = rlcout(:,1).*sin(slcout)

x_cone = [x_cone;xlc];
y_cone = [y_cone;ylc];


printf('mach: %f',mach(length(mach)))

psibuffl = atan(-(rlcout(:,2).*cos(slcout)-rlcout(:,1).*sin(slcout))./(rlcout(:,2).*sin(slcout)+rlcout(:,1).*cos(slcout)));
psibuffu =  atan((y_cone(1:19,1)-y_uc)./x_cone(1:19,1));
psi = [psibuffu;psibuffl]
figure(3)
for j = 1:rows(x_cone)
	%局所円錐衝撃波面の傾きψ
	
	if y_cone(j) >= y_rec
		R_el(j) = r_uc;
		
	else
		R_el(j) = rlcout(j-19,3);
	end
	%中心点
	cx(j) = fsolve(@(cx)circle_solve_x(cx,x_cone(j),y_cone(j),psi(j),R_el(j)),-10);
	cy(j) = tan(psi(j))*(cx(j)-x_cone(j))+y_cone(j);
	if cx(j) == x_cone(j)
		cy(j) = fsolve(@(cy)circle_solve_y(cy,x_cone(j),y_cone(j),psi(j),R_el(j)),-10);
		cx(j) = (cy(j)-y_cone(j))/tan(psi(j))+x_cone(j);
	end
	cz(j) = R_el(j)/sin(beta-delta);
	cycl(j) = tan(psi(j))*(0-x_cone(j))+y_cone(j);
	if j==1
		cycl(j)=y_uc;
	end
	czcl(j) = interp1([y_cone(j),cy(j)],[0,cz(j)],cycl(j),'extrap');
	plot([cx(j),x_cone(j)],[cy(j),y_cone(j)],'-o')
	hold on

end
axis 'equal'
grid on
hold off
