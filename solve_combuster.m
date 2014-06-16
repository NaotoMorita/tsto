function res = solve_combuster(f,Gio,Tio,Air_ent,H2_ent,Tmc,Aio,rhoio,Vio,T0io,pio,Ra)

	%水素燃料
	Rf = 4121.74 ;

	%主燃焼機空気流量
	Gamb = Gio;
	%空燃比
	Gfmb = Gamb/f ;
	
	
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
	qlf = 10.78/0.0899 * 1000000; %J/kg;
	Qfmb = etamb*qlf*Gfmbr
	Imbi = Gamb*interp1(Air_ent(:,1),Air_ent(:,2),Tio)*1000000 +  Gfmb*interp1(H2_ent(:,1),H2_ent(:,2),Tio)*1000000
	Imbo = Gamb*interp1(Air_ent(:,1),Air_ent(:,2),Tmc)*1000000 +  Gfmb*interp1(H2_ent(:,1),H2_ent(:,2),Tmc)*1000000
	
	Vio
	Rmc = (Ra *f+Rf)/(1+f);
	kappa_mc = hinetu_calc(f)
	rhomb2 = pio/Rmc/Tmc
	Vmb = (Gfmb+Gamb) / rhomb2 / Aio
	Gfmb
	Gamb
	res = Imbi+Qfmb + (Gfmb+Gamb) *Vio^2/2  -Imbo - (Gfmb+Gamb) *Vmb^2/2
	%asdf
end