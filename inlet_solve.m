function res = inlet_solve(MCR,TPR,pio_pin,Aio,Vin,Aii,rhoin,P0in,p0,T0in,kappa,Ra)

	%質量流量
	Gio = rhoin*Vin*Aii*MCR;
	%インテーク出口全圧
	P0io = P0in *TPR;



	%インテーク出口圧力
	pio =  p0*pio_pin;
	%インテーク出口マッハ数
	Mio = fsolve(@(M)P0M_solve(P0io,pio,M),2.5);
	%インテーク出口温度
	Tio = T0in/(1+(kappa-1)/2*Mio^2);
	%インテーク出口密度
	rhoio = pio/Ra/Tio;
	%インテーク出口速度
	Vio = Mio*sqrt(kappa*pio/rhoio);
	Vio2 = Mio*sqrt(kappa*Ra*Tio);
	Gio;
	Aiod = Gio/(Vio*rhoio);
	res =Aio-Aiod;
end