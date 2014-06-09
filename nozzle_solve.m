function res = nozzle_solve(ER,Mni,Mno,kappa_mc)
	res = Mni/Mno*((1+(kappa_mc-1)/2*Mno^2)/(1+(kappa_mc-1)/2*Mni^2))^((kappa_mc+1)/(2*(kappa_mc-1)))-ER;
end