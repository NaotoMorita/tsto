function [F,Isp,Gfmb] = scram_thrust(MCR,TPR,pio_pin,Tmc,Min,Tin,p0,rhoin,Aii,Aio,Ano)

	%---------------------------------------------------------インテーク


	kappa = 1.4;
	
	Ra = p0/Tin/rhoin;
	cin = sqrt(kappa*p0/rhoin);
	Vin = Min * cin;
	Dynamic_Pressure = 0.5*rhoin*Vin^2/1000
	P0in=p0*(1+(kappa-1)/2*Min^2)^((kappa)/(kappa-1));
	T0in = Tin*(1+(kappa-1)/2*Min^2);


	%質量流量
	Gio = rhoin*Vin*Aii*MCR;
	%インテーク出口全圧
	P0io = P0in *TPR;



	%インテーク出口圧力
	pio =  p0*pio_pin;
	%インテーク出口マッハ数
	Mio = fsolve(@(M)P0M_solve(P0io,pio,M),2.5)
	%インテーク出口温度
	Tio = T0in/(1+(kappa-1)/2*Mio^2);
	%インテーク出口密度
	rhoio = pio/Ra/Tio;
	%インテーク出口速度
	Vio = Mio*sqrt(kappa*pio/rhoio)
	Vio2 = Mio*sqrt(kappa*Ra*Tio)
	Gio
	Aio = Gio/(Vio*rhoio)

	%-----------------------------------------------------------燃焼機

		

		
	%非エンタルピー補間関数を作成
	H2_ent = [298 0;
			300 0.027;
			400 1.468;
			500 2.920;
			600 4.374;
			700 5.832;
			800 7.298;
			900 8.776;
			1000 10.27;
			1100 11.78;
			1200 13.30;
			1300 14.84;
			1400 16.41;
			1500 18.00;
			1600 19.62;
			1700 21.25;
			1800 22.91;
			1900 24.58;
			2000 26.27;
			2100 27.98;
			2200 29.71;
			2300 31.45;
			2400 33.21;
			2500 34.99;
			2600 36.78];
			
	Air_ent = [298 0;
			300 0.0019;
			400 0.1028;
			500 0.2051;
			600 0.3088;
			700 0.4151;
			800 0.5234;
			900 0.6347;
			1000 0.7478;
			1100 0.8637;
			1200 0.9797;
			1300 1.098;
			1400 1.216;
			1500 1.377;
			1600 1.459;
			1700 1.582;
			1800 1.705;
			1900 1.829;
			2000 1.952;
			2100 2.079;
			2200 2.205;
			2300 2.334;
			2400 2.458;
			2500 2.585;
			2600 2.714];
			
			

	%主燃焼機空気流量
	Gamb = Gio;


	
	%燃焼温度
	f = fsolve(@(f)solve_combuster(f,Gio,Tio,Air_ent,H2_ent,Tmc),71)
	if Tmc >= 2599.9
		printf("\n\nwarning!! Tmc over 2600K  outer of entalpy list \n\n\n")
	end
	Gfmb= Gamb/f
	%水素燃料
	Rf = 4121.74 ;



	%当量比
	phimb = 31.7/f;
	printf("phi = %.4f\n",phimb);
	%発熱燃料流量
	if phimb < 1
		Gfmbr = Gfmb;
	else 
		Gfmbr = Gfmb/phimb;
	end

	%主燃焼機発熱量
	etamb = 0.95;
	qlf = 10.78/0.0899; %MJ/m^3;
	Qfmb = etamb*qlf*Gfmbr;
	%エンタルピーにて水蒸気を考慮すること
	Imbi = Gamb*interp1(Air_ent(:,1),Air_ent(:,2),Tio)  +Gfmb*interp1(H2_ent(:,1),H2_ent(:,2),Tio);
	Imbo = Gamb*interp1(Air_ent(:,1),Air_ent(:,2),Tmc) +  Gfmb*interp1(H2_ent(:,1),H2_ent(:,2),Tmc)


	Rmc = (Ra *f+Rf)/(1+f);
	kappa_mc = hinetu_calc(f);

	%燃焼機音速
	cmb=sqrt(Tmc*kappa_mc*Rmc);
	%燃焼後密度
	rhomc = pio/Rmc/Tmc
	%燃焼機流量
	Gmb = Gamb+Gfmb;
	%燃焼後速度
	Vmc = Gmb/rhomc/Aio
	Vio
	Vio2 = Gmb/rhoio/Aio
	%燃焼後マッハ数
	Mmc = Vmc/cmb
	
	if Mmc<1
		printf("\n\nwarning!! Mmc:%.2f  ram conbustion\n\n\n",Mmc)
	end
	%全圧より静圧を求める
	pmc =  pio;
	P0mb=pmc*(1+(kappa_mc-1)/2*Mmc^2)^((kappa_mc)/(kappa_mc-1));
	T0mb = Tmc*(1+(kappa_mc-1)/2*Mmc^2);




	%排気ノズル
	%超音速燃焼のため、単純な等エントロピー膨張を用いてインテーク面積まで膨張させる
	Mni = Mmc;
	Ani = Aio;
	%Ano = Aii*0.93;
	ER = Ano/Ani;
	pni = pmc;
	Tni = Tmc;
	cni = cmb;

	Mno = fsolve(@(Mno)nozzle_solve(ER,Mni,Mno,kappa_mc),4);
	pno = pni * ((1+(kappa_mc-1)/2*Mni^2)/(1+(kappa_mc-1)/2*Mno^2))^(kappa_mc/(kappa_mc-1))
	Tno = Tni * ((1+(kappa_mc-1)/2*Mni^2)/(1+(kappa_mc-1)/2*Mno^2));
	cno = cni *  ((1+(kappa_mc-1)/2*Mni^2)/(1+(kappa_mc-1)/2*Mno^2))^(1/2);
	Vno = Mno * cno;

	%エンジン抗力係数をお求める
	Re = 1*Vin/1.5989*10^4;
	Cf = 0.472/(log10(Re)^2.58*(1+(kappa-1)/2*Min^2)^0.467);


	F = Gmb * (Vno-Vin) + (pno-p0)*Ano - 1/2*rhoin*Vin^2*Aii*Cf;
	printf("F(tonf)=%.3f\n",F/9.8/1000);
	Isp = F/9.8/Gfmb
	
	fid = fopen("Engine_Spec.txt","wt");
	fprintf(fid,"-----インテーク性能-----\n")
	fprintf(fid,"MCR : %.3f\n",MCR)
	fprintf(fid,"TPR : %.3f\n",TPR)
	fprintf(fid,"入口面積 : %.3f\n",Aii)	
	fprintf(fid,"流入マッハ数 : %.3f\n流入静圧 : %.3f\n流入静温 : %.3f\n流入密度 : %.3f\n流入空気ガス定数 : %.3f\n",Min,p0,Tin,rhoin,Ra)
	fprintf(fid,"流入総圧 : %.3f\n流入総温 : %.3f\n",P0in,T0in)
	fprintf(fid,"出口面積 : %.3f\n",Aio)
	fprintf(fid,"流出マッハ数 : %.3f\n流出静圧 : %.3f\n流出静温 : %.3f\n流出密度 : %.3f\n",Mio,pio,Tio,rhoio)
	fprintf(fid,"流出総圧 : %.3f\n流出総温 : %.3f\n",P0io,T0in)
	fprintf(fid,"-----燃焼機性能-----\n")
	fprintf(fid,"空気流量 : %.3f\n",Gamb)
	fprintf(fid,"燃料流量 : %.3f\n",Gfmb)
	fprintf(fid," 空燃比: %.3f\n",f)
	fprintf(fid,"当量比 : %.3f\n",phimb)
	fprintf(fid,"燃焼温度 : %.3f\n",Tmc)
	fprintf(fid,"燃焼後比熱 : %.3f\n",kappa_mc)
	fprintf(fid,"流出マッハ数 : %.3f\n流出静圧 : %.3f\n流出静温 : %.3f\n流出ガス定数 : %.3f\n",Mmc,pmc,Tmc,Rmc)
	fprintf(fid,"燃焼後総圧 : %.3f\n",P0mb)
	fprintf(fid,"燃焼後全温 : %.3f\n",T0mb)
	fprintf(fid,"-----ノズル性能-----\n")
	fprintf(fid,"ノズル入口面積 : %.3f\n",Ani)
	fprintf(fid,"ノズル開口比: %.3f\n",ER)
	fprintf(fid,"ノズル出口面積 : %.3f\n",Ano)
	fprintf(fid,"流出マッハ数 : %.3f\n流出静圧 : %.3f\n流出静温 : %.3f\n",Mno,pno,Tno)
	fprintf(fid,"機体抗力 : %.3f\n",1/2*rhoin*Vin^2*Aii*Cf)
	fprintf(fid,"推力 : %.3f\n",F)
	fprintf(fid,"比推力 : %.3f\n",Isp)
	fclose(fid);
	
end