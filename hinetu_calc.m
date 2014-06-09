function hinetu = hinetu_calc(f)
	if f <= 14
		hinetu = 2.813e-6*f^5-1.16e-4*f^4+3.598e-3*f^3-3.831e-2*f^2+1.859e-1*f+1.009e0;
	else
		hinetu = 4.492e-14*f^6-4.232e-11*f^5+1.560e-8*f^4-2.830e-6*f^3+2.585e-4*f^2-1.005e-2*f+1.376e0 ;
	end
end