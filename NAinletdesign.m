1;
clear all
function res = prandtl_meyer(theta, mach, beta, kappa)
	res = tan(theta) - ((mach ^ 2 * sin(beta) ^ 2 - 1) * cot(beta) / (1 + (1 / 2 * (kappa + 1) - sin(beta) ^ 2) * mach^2));
endfunction

function res = PMsolve(gamma, mach, beta, kappa)
	res = tan(beta-gamma) - ((mach ^ 2 * sin(beta) ^ 2 - 1) * cot(beta) / (1 + (1 / 2 * (kappa + 1) - sin(beta) ^ 2) * mach^2));
endfunction

function mach2 = PMmach(mach, beta, kappa)
	ms = mach * sin(beta);
	mach2     = ((1.0 + (kappa - 1.0) * ms ^2 + mach ^ 4.0 * (1.0 / 4.0 * (kappa + 1.0) ^ 2.0 - kappa * sin(beta) ^ 2.0) * sin(beta) ^ 2.0) / ((1.0 + 1.0 / 2.0 * (kappa - 1.0) * ms ^ 2.0) * (kappa * ms ^ 2.0 - 1.0 / 2.0 * (kappa - 1.0)))) ^ (1.0 / 2.0);
end


%主流マッハ数
M0 = 6
%最終φ > 90 deg
phiend = deg2rad(85);
%初期φ　
phi0 = deg2rad(20);
kappa = 1.4;

dphi = 0.0001

%初期theta
theta0 = fsolve(@(theta)prandtl_meyer(theta,M0,phi0,kappa),0.0001);
theta(1) = theta0;
mach(1) = PMmach(M0,phi0,kappa);
beta(1) = phi0;
T(1) = theta0;
phi(1) = phi0;

for i = 2:100
	phi(i) = phi(i-1) + dphi;
	beta(i) = phi(i)-theta(i-1);
	T(i) =  fsolve(@(theta)prandtl_meyer(theta,mach(i-1),beta(i),kappa),0.0001);
	theta(i) = T(i) + theta(i-1);
	mach(i) =  PMmach(mach(i-1),beta(i),kappa);
end
