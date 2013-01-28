clc;
clear('all');

Nden = 2E3;

% Triple line generation
Tt = TriplePointTemperature();
[rhotl,rhotg] = TriplePointDensities();
rhot = GeometricSpace(rhotg,rhotl,1.02,Nden);
it = InternalEnergy(rhot,Tt);
Pt = Pressure(rhot,Tt);

% Saturation line generation
[Psat,Tsat,rhol,rhog] = SaturationStateGivenDensity(rhot);
rhoSat = [rhog;rhol(end:-1:1)];
iSat = InternalEnergy(rhoSat,[Tsat;Tsat(end:-1:1)]);

% Isotherm generation
Tisos  = [300,373.15,393,450,500,550,647.096];
[rhoIso,Tisos] = meshgrid(rhot,Tisos);
IntIso = InternalEnergy(rhoIso,Tisos);
Piso   = Pressure(rhoIso,Tisos);

Mask   = [false,true,true,false,false,false,false];

figure(1);
semilogx(1./rhot,it/1E6,'y','LineWidth',2);
hold('on');
semilogx(1./rhoSat,iSat/1E6,'y','LineWidth',2);
for k = 1:size(Tisos,1)
    if not(Mask(k))
        semilogx(1./rhot,IntIso(k,:)'/1E6,'b','LineWidth',2);
    else
        semilogx(1./rhot,IntIso(k,:)'/1E6,'b--','LineWidth',2);
    end
end
axis([5E-4,3E2,-0.25,3.2]);
grid('on');
xlabel('Volume [m^3/kg]');
ylabel('Internal Energy [MJ/kg]');
title('Wisconsin Water Property Package');
hold('off');
MakePNG('InternalEnergyVsVolume');

figure(2);
Patm = 101325;
loglog(1./rhot,Pt/Patm,'k','LineWidth',2);
hold('on');
loglog(1./rhoSat,[Psat;Psat(end:-1:1)]/Patm,'k','LineWidth',2);
for k = 1:size(Tisos,1)
    if not(Mask(k))
        loglog(1./rhot,Piso(k,:)'/Patm,'b','LineWidth',2);
    else
        loglog(1./rhot,Piso(k,:)'/Patm,'b--','LineWidth',2);
    end
end
axis([5E-4,3E2,2E-3,3E3]);
xlabel('Volume [m^3/kg]');
ylabel('Pressure [atm]');
grid('on');
title('Wisconsin Water Property Package');
hold('off');
MakePNG('PressureVsVolume');

