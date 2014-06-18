
%OCBIを形作る
	for j = 1:rows(x_cone)
		if x_cone(j)>x_rec
			xstt(j) = x_rec;
			ystt(j) = tan(psi(j))*(x_rec-x_cone(j))+y_cone(j);
			r_rec(j) = sqrt((xstt(j)-cx(j))^2+(ystt(j)-cy(j))^2);
		else
			%ystt(j) = y_cone(j);
			%xstt(j) = x_cone(j);
			ystt(j) = y_rec;
			xstt(j) = (y_rec-y_cone(j))/tan(psi(j))+x_cone(j);
			r_rec(j) = sqrt((xstt(j)-cx(j))^2+(ystt(j)-cy(j))^2);
		end
		zstt(j) = (R_el(j)-r_rec(j))/tan(beta-delta);
		rcone(j) =r_rec(j)/sin(beta-delta);
		x_oc(:,j) = xstt(j)+rcone(j).*(x'.*cos(psi(j))); 
		y_oc(:,j) = ystt(j)+rcone(j).*(x'.*sin(psi(j))); 
		z_oc(:,j) = -zstt(j)+rcone(j).*R_el(j).*(z');
	end


figure(1)
nowwork=pwd;
cd OCBI
for j = 1:rows(x_cone)
	plot3(x_oc(:,j),y_oc(:,j),z_oc(:,j))
	hold on
	if cx(j)<0.0001
		plot3([xstt(j),0],[ystt(j),cycl(j)],-[zstt(j),czcl(j)],'-r')
		plot3(-[xstt(j),0],[ystt(j),cycl(j)],-[zstt(j),czcl(j)],'-r')
	else
		plot3([xstt(j),cx(j)],[ystt(j),cy(j)],-[zstt(j),cz(j)],'-r')
		plot3(-[xstt(j),cx(j)],[ystt(j),cy(j)],-[zstt(j),cz(j)],'-r')
	end
	z_out = linspace(z_oc(1,j),z_oc(rows(x_oc),j),200)';
	x_out = spline(z_oc(:,j),x_oc(:,j),z_out);
	y_out = spline(z_oc(:,j),y_oc(:,j),z_out);
	if min(x_oc(:,j)) <= 0.0001
		dlmwrite(strcat('stream',mat2str(j),'.sldcrv'),[zeros(rows(x_out),1),y_out,z_out].*1000);
	elseif min(y_oc(:,j)) <= 0.0001
		dlmwrite(strcat('stream',mat2str(j),'.sldcrv'),[x_out,zeros(rows(y_out),1),z_out].*1000);
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
for j = 1:rows(x_cone)
	if cx(j)<0.0001
		plot3([0,xstt(j)],[cycl(j),ystt(j)],[czcl(j),zstt(j)],'r-o')
		hold on
		plot3(-[0,xstt(j)],[cycl(j),ystt(j)],[czcl(j),zstt(j)],'r-o')
	else
		plot3([cx(j),xstt(j)],[cy(j),ystt(j)],[cz(j),zstt(j)],'r-o')
		hold on
		plot3(-[cx(j),xstt(j)],[cy(j),ystt(j)],[cz(j),zstt(j)],'r-o')
	end
	plot3([x_cone(j),xstt(j)],[y_cone(j),ystt(j)],[0,zstt(j)])
	plot3(-[x_cone(j),xstt(j)],[y_cone(j),ystt(j)],[0,zstt(j)])
	
end
axis 'equal'
view(210,45)
grid on
hold off