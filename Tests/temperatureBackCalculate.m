% clc();
clear();


%   Quality and tau
Nx  = 50;
Nt  = 50;
x   = repmat(linspace(0,1,Nx)',Nt,1);
tau = linspace(1.000001*CriticalTemperatureR(),TriplePointTau(),Nt)';

%   Saturation state
[~,delL,delG] = SaturationStateGivenTauRRND(tau);
iL = InternalEnergyOneRND(delL,tau);
iG = InternalEnergyOneRND(delG,tau);

%   Calculate mixutre energies
I    = ceil((1:(Nx*Nt))'/Nx)                ;
iMix = x.*iL(I) + (1-x).*iG(I)              ;
dMix = 1./(x./delL(I) + (1-x).*1./delG(I))  ;

%   Back calculate taus
tic;
[~,tauB,~,~,~] = SaturationStateGivenMixedRhoIRRND(dMix,iMix);
toc;

%   Plot
figure(1);
    semilogy((1:numel(tauB))',abs(tau(I) - tauB));

