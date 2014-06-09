1;
clear all
function res = P0M_solve(P0,p,M)
	kappa = 1.4;
	res = p*(1+(kappa-1)/2*M^2)^((kappa)/(kappa-1))-P0;
end

%---------------------------------------------------------インテーク
TPR = 0.60
MCR = 1.4;
%空気
Ra = 287.04;

%主流マッハ数
Min = 6

%主流音速
kappa = 1.4

Tin = 216
rhoin=0.088
p0 = 5531
Ra = p0/Tin/rhoin
cin = sqrt(kappa*p0/rhoin)
Vin = Min * cin;
P0in=p0*(1+(kappa-1)/2*Min^2)^((kappa)/(kappa-1))
T0in = Tin*(1+(kappa-1)/2*Min^2)


%インテーク取り込み面積
Aii = 3.66*2.02;
%インテークスロート面積
Aio = 0.28*2.02;

%質量流量
Gio = rhoin*Vin*Aii*MCR;
%インテーク出口全圧
P0io = P0in *TPR

%インテーク圧縮比
pio_pin = 50;
%インテーク出口速度
Vio = 4*(P0io-pio_pin*p0)*Aio/Gio
%インテーク出口密度
rhoio = Gio/Vio/Aio
%インテーク出口温度


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





%-----------------------------------------------------------燃焼機

function Tmc =Tmc_calc(f)
	if f <= 50
		Tmc = 3.612e-5*f^5-4.566e-3*f^4+2.209e-1*f^3-6.874e0*f^2+1.691e2*f+3.155e2;
	else
		Tmc = -7.569e-9*f^5+7.308e-6*f^4-2.814e-3*f^3+5.543e-1*f^2-5.983e1*f+3.855e3 ;
	end
end
	
function hinetu = hinetu_calc(f)
	if f <= 14
		hinetu = 2.813e-6*f^5-1.16e-4*f^4+3.598e-3*f^3-3.831e-2*f^2+1.859e-1*f+1.009e0;
	else
		hinetu = 4.492e-14*f^6-4.232e-11*f^5+1.560e-8*f^4-2.830e-6*f^3+2.585e-4*f^2-1.005e-2*f+1.376e0 ;
	end
end
	
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
		
		

%空燃比
f = 35;
function res = solve_combuster(f,Gio,Tio,Air_ent,H2_ent,Tmc)

	%燃焼温度
	%Tmc = Tmc_calc(f);
	%水素燃料
	Rf = 4121.74 ;

	%主燃焼機空気流量
	Gamb = Gio;
	%主燃焼機空気流量
	Gfmb = Gamb/f;

	%当量比
	phimb = 31.7/f;
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
	Imbi = Gamb*interp1(Air_ent(:,1),Air_ent(:,2),Tio) +  Gfmb*interp1(H2_ent(:,1),H2_ent(:,2),Tio);
	Imbo = Gamb*interp1(Air_ent(:,1),Air_ent(:,2),Tmc) +  Gfmb*interp1(H2_ent(:,1),H2_ent(:,2),Tmc);

	res = Imbi+Qfmb-Imbo;
	
end
%燃焼温度
f = 60%Tmc_calc(f)
Tmc = fsolve(@(Tmc)solve_combuster(f,Gio,Tio,Air_ent,H2_ent,Tmc),2000);

%水素燃料
Rf = 4121.74 ;

%主燃焼機空気流量
Gamb = Gio;
%主燃焼機空気流量
Gfmb = Gamb/f;

%当量比
phimb = 31.7/f;
%発熱燃料流量
if phimb < 1
	Gfmbr = Gfmb;
else 
	Gfmbr = Gfmb/phimb;
end

%主燃焼機発熱量

etamb = 0.95;
qlf = 10.78/0.0899; %MJ/m^3;
Qfmb = etamb*qlf*Gfmbr
Imbi = Gamb*interp1(Air_ent(:,1),Air_ent(:,2),Tio) +  Gfmb*interp1(H2_ent(:,1),H2_ent(:,2),Tio)
Imbo = Gamb*interp1(Air_ent(:,1),Air_ent(:,2),Tmc) +  Gfmb*interp1(H2_ent(:,1),H2_ent(:,2),Tmc)


Rmc = (Ra *f+Rf)/(1+f);
kappa_mc = hinetu_calc(f);

%燃焼機音速
cmb=sqrt(Tmc*kappa_mc*Rmc)
%速度は変化無しだから
Mmc = Vio/cmb;

%全圧より静圧を求める
pmc =  pio;
P0mb=pmc*(1+(kappa_mc-1)/2*Mmc^2)^((kappa_mc)/(kappa_mc-1))
T0mb = Tmc*(1+(kappa_mc-1)/2*Mmc^2)

%P0mb 主燃焼機全圧
%Tmc 主燃焼機燃焼静温
%kappa_mc 主燃焼機比熱比

%燃焼機流量
Gmb = Gamb+Gfmb;




%排気ノズル
%超音速燃焼のため、単純な等エントロピー膨張を用いてインテーク面積まで膨張させる
Mni = Mmc;
Ani = Aio;
Ano = Aii*0.93
ER = Ano/Ani;
pni = pmc;
Tni = Tmc;
cni = cmb;

function res = nozzle_solve(ER,Mni,Mno,kappa_mc)
	res = Mni/Mno*((1+(kappa_mc-1)/2*Mno^2)/(1+(kappa_mc-1)/2*Mni^2))^((kappa_mc+1)/(2*(kappa_mc-1)))-ER;
end

Mno = fsolve(@(Mno)nozzle_solve(ER,Mni,Mno,kappa_mc),4)
pno = pni * ((1+(kappa_mc-1)/2*Mni^2)/(1+(kappa_mc-1)/2*Mno^2))^(kappa_mc/(kappa_mc-1))
Tno = Tni * ((1+(kappa_mc-1)/2*Mni^2)/(1+(kappa_mc-1)/2*Mno^2))
cno = cni *  ((1+(kappa_mc-1)/2*Mni^2)/(1+(kappa_mc-1)/2*Mno^2))^(1/2)
Vno = Mno * cno

%エンジン抗力係数をお求める
Re = 15*Vin/1.5989*10^4
Cf = 0.472/(log10(Re)^2.58*(1+(kappa-1)/2*Min^2)^0.467)


F = Gmb * (Vno-Vin) + (pno-p0)*Ano - 1/2*rhoin*Vin^2*Aii*Cf

%~ %臨界圧力
%~ pncr = P0mb*(2/(kappa_mc+1))^(kappa_mc/(kappa_mc-1));

%~ %スロート静圧
%~ if p0 >=pncr
	%~ pnt = p0;
%~ else
	%~ pnt = pncr;
%~ end

%~ %スロート静温
%~ Tnt = Tmc*(pnt/P0mb)^((kappa_mc-1)/kappa_mc);
%~ %スロート密度
%~ rhont = pnt/Rmc/Tnt;

%~ %ノズルスロート速度
%~ Cnt = sqrt(kappa_mc*Rmc*Tnt);
%~ if p0 >= pncr
	%~ Vnt = sqrt(2*Rmc*Tmc*kappa_mc/(kappa_mc-1)*(1-(p0/P0mb)^((kappa_mc-1)/kappa_mc)));
%~ else
	%~ Vnt = Cnt; 
%~ end

%~ %ノズルスロート流量
%~ Gnt = Gmb; 
%~ Ant = Gnt/rhont/Vnt;

%~ %ノズルスロートマッハ数
%~ Mnt = Vnt/Cnt;

%~ %全温の保存より膨張比を求める
%~ function res = epn_solve(epn,Tnt,Mnt,kappa_mc,T0mb)
	%~ Tno = Tnt*(1/epn)^((kappa_mc-1)/kappa_mc)
	%~ if Mnt  < 1
		%~ Mno = Mnt;
		%~ ER=1
	%~ else
		%~ Mno = sqrt(2/(kappa_mc-1) * epn ^ ((kappa_mc-1)/kappa_mc))
		%~ ER = 1/Mno*(((kappa_mc-1)*Mno^2+2)/(kappa_mc-1))^((kappa_mc+1/(2*(kappa_mc-1))));
	%~ end
	%~ T0no = Tno*(1+(kappa_mc-1)/2*Mno^2)

	%~ res = T0no-T0mb
%~ end
	
%~ epn = fsolve(@(epn)epn_solve(epn,Tnt,Mnt,kappa_mc,T0mb),0.1);
%~ if Mnt  < 1
	%~ Mno = Mnt;
	%~ ER=1
%~ else
	%~ Mno = sqrt(2/(kappa_mc-1) * epn ^ ((kappa_mc-1)/kappa_mc))
	%~ ER = 1/Mno*(((kappa_mc-1)*Mno^2+2)/(kappa_mc-1))^((kappa_mc+1/(2*(kappa_mc-1))));
%~ end

%~ function Mno_solve

%~ if ER > Aii/Ant;
	
%~ Ano = Ant*ER