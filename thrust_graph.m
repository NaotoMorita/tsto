1;
for i = 1:11
	TPR(i) = 0.5+0.05*(i-1);
	for j = 1:11
		p_pin(j) = 10+4*(j-1);
		Aio(j) = (-(1.14-0.28)/10*(j-1)+1.14)*2.02;
		[F(i,j) Isp(i,j) Gfmb(i,j)] = scram_thrust(1,TPR(i),p_pin(j),1900,6,216,5531,0.088,3.66*2.02,Aio(j),3.66*2.02);
	end
end
	
figure(1)
	mesh(TPR,p_pin,F/1000)
	ax=xlabel("TPR")
	set(ax,'fontsize',15)
	ay=ylabel("p/pin")
	set(ay,'fontsize',15)
	az=zlabel("F[kN]")
	set(az,'fontsize',15)
	print "F.png" -dpng
	
figure(2)
	mesh(TPR,p_pin,Isp)
	ax=xlabel("TPR")
	set(ax,'fontsize',15)
	ay=ylabel("p/pin")
	set(ay,'fontsize',15)
	az=zlabel("Isp")
	set(az,'fontsize',15)
	print "Isp.png" -dpng
	
figure(3)
	mesh(TPR,p_pin,Gfmb)
	ax=xlabel("TPR")
	set(ax,'fontsize',15)
	ay=ylabel("p/pin")
	set(ay,'fontsize',15)
	az=zlabel("Gfmb")
	set(az,'fontsize',15)
	print "Gfmb.png" -dpng