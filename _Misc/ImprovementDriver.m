function [] = ImprovementDriver(rho,Scale)
    
[~,Tsat,~,~] = SaturationStateGivenDensity(rho);


% Show([Psat,tauSat,delL,delG]);
% [delLt,delGt] = TriplePointDensities();
Tc   = CriticalTemperature();
Tt   = TriplePointTemperature();


N = 500;
Tspace = linspace(Tt,Tsat,N);
IntEnergy = InternalEnergy(rho,Tspace);


it   = IntEnergy(1);
isat = IntEnergy(end);

IndexMid = round(Scale*N);
Tmid = Tspace(IndexMid);
imid = IntEnergy(IndexMid);


iFit = LagrangeInterpolation(linspace(0,1,N),[0,Scale,1],[it,imid,isat]);

figure();
plot(Tc./Tspace,IntEnergy/1E6,Tc./Tspace,iFit/1E6,'r--',Tc./[Tt,Tmid,Tsat],[it,imid,isat]/1E6,'g.');
title(['\rho = ',num2str(rho)]);

    
    
end