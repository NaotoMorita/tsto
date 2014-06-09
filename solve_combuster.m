function res = solve_combuster(f,Gio,Tio,Air_ent,H2_ent,Tmc)

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
	qlf = 10.78/0.0899; %MJ/m^3;
	Qfmb = etamb*qlf*Gfmbr;
	Imbi = Gamb*interp1(Air_ent(:,1),Air_ent(:,2),Tio) +  Gfmb*interp1(H2_ent(:,1),H2_ent(:,2),Tio);
	Imbo = Gamb*interp1(Air_ent(:,1),Air_ent(:,2),Tmc) +  Gfmb*interp1(H2_ent(:,1),H2_ent(:,2),Tmc);

	res = Imbi+Qfmb-Imbo;
	
end