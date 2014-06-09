%流れから、ショック角と転向角、マッハ数
1;
clear all
function res = prandtl_meyer(theta, mach, beta, kappa)
	res = tan(theta) - ((mach ^ 2 * sin(beta) ^ 2 - 1) * cot(beta) / (1 + (1 / 2 * (kappa + 1) - sin(beta) ^ 2) * mach^2));
endfunction
	
%main関数
kappa = 1.4;

want = menu("waht to want","mach angle","amplitude angle","mach number");

if want == 1
	theta = deg2rad(input("input amplitude angle(deg) "));
	mach = input("input mach ");
	beta = rad2deg(fsolve(@(beta)prandtl_meyer(theta,mach,beta,kappa),0.0001))
elseif want == 2
	mach = input("input mach ");
	beta = deg2rad(input("input mach angle(deg) "));
	theta = rad2deg(fsolve(@(theta)prandtl_meyer(theta,mach,beta,kappa),0.0001))
else 
	beta = deg2rad(input("input mach angle(deg) "));
	theta = deg2rad(input("input amplitude angle(deg) "));
	mach = fsolve(@(mach)prandtl_meyer(theta,mach,beta,kappa),0.0001)
end