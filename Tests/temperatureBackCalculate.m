% clc();
clear();



% [~,delL]   = saturableDeltas();
% [~,tau]    = tausAtSaturableDeltas();
tau           = 2.367;
[~,delL,delG] = SaturationStateGivenTauRRND(tau);
delta      = logspace(log10(delG),log10(delL),10E3)';
tau        = delta*0 + tau;
iND        = InternalEnergy(delta*322,647.096./tau)/DimensioningInternalEnergy();

tic;
[~,tauB,~,~,~] = SaturationStateGivenMixedRhoIRRND(delta,iND);
toc;
any(isnan(tauB))
semilogy((1:10E3)',abs(tau - tauB));

