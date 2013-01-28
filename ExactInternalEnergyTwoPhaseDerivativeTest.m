clc;
clear('all');

T = 300;
tau = CriticalTemperature./T;
[Psat,rhol,rhog] = SaturationStateGivenTsat(T);

rho = 0.9999*rhol;
del  = rho  / CriticalDensity();
delL = rhol / CriticalDensity();
delG = rhog / CriticalDensity();

dT = 1E-5;
i_TN = (InternalEnergy(rho,T+dT)-InternalEnergy(rho,T-dT))./(2*dT);
tau_T = tau_Temperature(T);
i_tauN = i_TN./tau_T;

il      = InternalEnergyOne(rhol,T);
ig      = InternalEnergyOne(rhog,T);
il_del  = InternalEnergyOne_delta(rhol,T);
ig_del  = InternalEnergyOne_delta(rhog,T);
il_tau  = InternalEnergyOne_tau(rhol,T);
ig_tau  = InternalEnergyOne_tau(rhog,T);

x       = QualityFromDensity(rho,rhol,rhog);
Pnd     = Psat/DimensioningPressure();
Pnd_tau = ClausiusClapeyronRRND(Psat,tau,delL,delG,il,ig);
PhiR_d  = HelmholtzResidual_d ([delL;delG],tau);
PhiR_dd = HelmholtzResidual_dd([delL;delG],tau);
PhiR_dt = HelmholtzResidual_dt([delL;delG],tau);
dl_tau  = (Pnd + tau.*Pnd_tau-delL.^2.*PhiR_dt(1))./(1+2.*delL.*PhiR_d(1)+delL.^2.*PhiR_dd(1));
dg_tau  = (Pnd + tau.*Pnd_tau-delG.^2.*PhiR_dt(2))./(1+2.*delG.*PhiR_d(2)+delG.^2.*PhiR_dd(2));

x_tau = ((delL-del).*delL.*dg_tau + (del-delG).*delG.*dl_tau)./(del.*(delG-delL).^2);

i_tauExact = x_tau.*(il - ig) + (1-x).*(ig_tau + dg_tau.*ig_del) + x.*(il_tau + dl_tau.*il_del);
i_TExact   = -i_tauExact.*tau_T;

Show(i_tauN - i_tauExact);
Show((i_tauN - i_tauExact)./i_tauExact);


subplot(2,1,1);
Tplot = linspace(-10,10)+T;
iPlot = InternalEnergy(rho,Tplot);
iTangent = i_TN*(Tplot-T)+InternalEnergy(rho,T);
plot(Tplot,iPlot,Tplot,iTangent,'ro');

subplot(2,1,2);
Tplot = linspace(-10,10)+T;
iPlot = InternalEnergy(rho,Tplot);
iTangent = i_TExact*(Tplot-T)+InternalEnergy(rho,T);
plot(Tplot,iPlot,Tplot,iTangent,'ro');

