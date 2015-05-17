clc();
clear();



[~,delL]   = saturableDeltas();
[~,tau]    = tausAtSaturableDeltas();
[~,~,delG] = SaturationStateGivenTauRRND(tau);
delta      = logspace(log10(delG),log10(delL),10E3)';
tau        = delta*0 + tau;
iND        = InternalEnergy(delta*322,647.096./tau)/DimensioningInternalEnergy();

[~,tauB,~,~,~] = SaturationStateGivenMixedRhoIRRND(delta,iND);