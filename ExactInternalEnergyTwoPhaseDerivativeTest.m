% clc;
clear('all');

T = 300;
tau = CriticalTemperature./T;
[Psat,rhol,rhog] = SaturationStateGivenTsat(T);

rho = 100 * rhog;
del  = rho  / CriticalDensity();
delL = rhol / CriticalDensity();
delG = rhog / CriticalDensity();

dT = 1E-9;
i_TN = (InternalEnergy(rho,T+dT)-InternalEnergy(rho,T-dT))./(2*dT);
tau_T = tau_Temperature(T);
i_tauN = i_TN./tau_T;

iL      = InternalEnergyOne(rhol,T);
iG      = InternalEnergyOne(rhog,T);

x        = QualityFromDensity(rho,rhol,rhog);
Pnd      = Psat/DimensioningPressure();
Pnd_tau  = ClausiusClapeyronRRND(Psat,tau,delL,delG,iL,iG);
PhiR_d   = HelmholtzResidual_d ([delL;delG],tau);
PhiR_dd  = HelmholtzResidual_dd([delL;delG],tau);
PhiR_dt  = HelmholtzResidual_dt([delL;delG],tau);
delL_tau = (delL + tau.^2 .* Pnd_tau + del.^2 .* (PhiR_d(1) - tau .* PhiR_dt(1))) ./ ...
           (tau .* (1 + 2 * del .* PhiR_d(1) + del.^2 .* PhiR_dd(1)));
delG_tau = (delG + tau.^2 .* Pnd_tau + del.^2 .* (PhiR_d(2) - tau .* PhiR_dt(2))) ./ ...
           (tau .* (1 + 2 * del .* PhiR_d(2) + del.^2 .* PhiR_dd(2)));
       
iL_del  = InternalEnergyOne_delta(rhol,T);
iG_del  = InternalEnergyOne_delta(rhog,T);
iL_tau  = iL_del .* delL_tau + InternalEnergyOne_tau(rhol,T);
iG_tau  = iG_del .* delG_tau + InternalEnergyOne_tau(rhog,T);

x_tau = ((delL-del).*delL.*delG_tau + (del-delG).*delG.*delL_tau)./(del.*(delG-delL).^2);

i_tauExact = iL_tau + x .* (iG_tau - iL_tau) + x_tau .* (iG - iL);
i_TExact   = i_tauExact.*tau_T;

drho = 1E-9;
i_rhoN  = (InternalEnergy(rho + drho,T)-InternalEnergy(rho - drho,T))./(2*drho);
rho_del =  CriticalDensity();
i_delN  = i_rhoN .* rho_del;

x_del = delG .* delL / (del.^2 .* (delG - delL));
i_del = x_del .* (iG - iL);
i_rho = i_del / CriticalDensity();


Show([(i_tauN - i_tauExact);(i_tauN - i_tauExact)./i_tauExact]);
Show([(i_delN - i_del);(i_delN - i_del)/i_del]);


subplot(2,1,1);
Tcen  = T;
iCen  = InternalEnergy(rho,Tcen);
Tplot = (0.999999:0.00000001:1.000001)' * Tcen;
iPlot = InternalEnergy(rho,Tplot);
plot(Tplot,iPlot,Tcen,iCen,'o',Tplot,i_TExact*(Tplot-Tcen) + iCen);


subplot(2,1,2);
rhoCen = rho;
iCen  = InternalEnergy(rhoCen,T);
rhoPlot = (0.99:0.01:1.01)' * rho;
iPlot = InternalEnergy(rhoPlot,T);
plot(rhoPlot,iPlot,rhoCen,iCen,'o',rhoPlot,i_rho*(rhoPlot-rhoCen) + iCen);



