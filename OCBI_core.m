%OCBI設計プログラム
1;
clear all

%パラメータ設定
%比熱比
kappa = 1.4
%衝撃波角
beta = deg2rad(34.3797988108415)
%衝撃波前マッハ数
M1 = 2.5
delta = atan(2*cot(beta)*((M1^2*sin(beta)^2-1)/(M1^2*(kappa+cos(2*beta))+2)));
Mn1 = M1*sin(beta);
Mn2 = sqrt((Mn1^2+(2/(kappa-1)))/(2*kappa*Mn1^2/(kappa-1)-1));
M2 = Mn2/sin(beta-delta);


%function読込
OCBI_functions

%Taylor-maccollを解く
OCBI_TM

%衝撃波設計
OCBI_shock

%流線とplot
OCBI_streamline

