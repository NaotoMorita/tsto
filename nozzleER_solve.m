function res = nozzleER_solve(Ano,Ani,Mni,pmc,p0,Tmc,cmb,kappa_mc)
	ER = Ano/Ani;
	pni = pmc;
	Tni = Tmc; 
	cni = cmb;

	Mno = fsolve(@(Mno)nozzle_solve(ER,Mni,Mno,kappa_mc),4);
	pno = pni * ((1+(kappa_mc-1)/2*Mni^2)/(1+(kappa_mc-1)/2*Mno^2))^(kappa_mc/(kappa_mc-1));
	res =pno-p0;
end