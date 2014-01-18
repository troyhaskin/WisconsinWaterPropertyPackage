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
i_tauN = i_TN./tau_T/DimensioningInternalEnergy();

iL      = InternalEnergyOneRND(delL,tau);
iG      = InternalEnergyOneRND(delG,tau);

x        = QualityFromDensity(rho,rhol,rhog);
Pnd      = Psat/DimensioningPressure();
Pnd_tau  = ClausiusClapeyronRRND(Pnd,tau,delL,delG,iL,iG);
PhiR_d   = HelmholtzResidual_d ([delL;delG],tau);
PhiR_dd  = HelmholtzResidual_dd([delL;delG],tau);
PhiR_dt  = HelmholtzResidual_dt([delL;delG],tau);
delL_tau = (delL + tau.^2 .* Pnd_tau + delL.^2 .* (PhiR_d(1) - tau .* PhiR_dt(1))) ./ ...
           (tau .* (1 + 2 * delL .* PhiR_d(1) + delL.^2 .* PhiR_dd(1)));
delG_tau = (delG + tau.^2 .* Pnd_tau + delG.^2 .* (PhiR_d(2) - tau .* PhiR_dt(2))) ./ ...
           (tau .* (1 + 2 * delG .* PhiR_d(2) + delG.^2 .* PhiR_dd(2)));
       
iL_del  = InternalEnergyOneRND_delta(delL,tau);
iG_del  = InternalEnergyOneRND_delta(delG,tau);
iL_tau  = iL_del .* delL_tau + InternalEnergyOneRND_tau(delL,tau);
iG_tau  = iG_del .* delG_tau + InternalEnergyOneRND_tau(delG,tau);

x_tau = ((delL-del).*delL.*delG_tau + (del-delG).*delG.*delL_tau)./(del.*(delG-delL).^2);

i_tauExact = iL_tau + x .* (iG_tau - iL_tau) + x_tau .* (iG - iL);
i_TExact   = i_tauExact.*tau_T;

drho = 1E-9;
i_rhoN  = (InternalEnergy(rho + drho,T)-InternalEnergy(rho - drho,T))./(2*drho);
rho_del =  CriticalDensity();
i_delN  = i_rhoN .* rho_del/DimensioningInternalEnergy();

x_del = delG .* delL / (del.^2 .* (delG - delL));
i_del = x_del .* (iG - iL);
i_rho = i_del / CriticalDensity();


Show([(i_tauN - i_tauExact);(i_tauN - i_tauExact)./i_tauExact]);
Show([(i_delN - i_del);(i_delN - i_del)/i_del]);




