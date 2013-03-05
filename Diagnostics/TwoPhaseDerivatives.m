clc;
clear('all');

delta = 10.^((-16:0)');
N     = length(delta);

rho  = 550   ;
T    = linspace(275,640,500);
Ints = InternalEnergy(rho,T);
dPdi = Pressure_InternalEnergy(rho,T);

[Psat,Tsat,rhoL,rhoG] = SaturationStateGivenDensity(rho);
isat    = InternalEnergy(rho,Tsat);
dPdiSat = Pressure_InternalEnergy(rho,Tsat);

figure(1);
title('P_i vs i: Absolute, Central Finite Differencing');
semilogx(Ints,dPdi,isat,dPdiSat,'ro');

