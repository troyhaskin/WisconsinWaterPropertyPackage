% clc();
clear();


%   Quality and tau
Nx  = 200;
Nt  = 600;
[delLt,delGt] = TriplePointDensitiesR();
% dMix = 10.^GeometricSpace(log10(delGt),log10(delLt),1.17,Nx);
dMix = GeometricSpace(delGt,delLt,1.13,Nx);
x    = repmat(QualityFromDensity(dMix,delLt,delGt),Nt,1);
tau  = linspace(1.00001*CriticalTemperatureR(),0.99999*TriplePointTau(),Nt)';

%   Saturation state
[~,delL,delG] = SaturationStateGivenTauRRND(tau);
iL = InternalEnergyOneRND(delL,tau);
iG = InternalEnergyOneRND(delG,tau);

%   Calculate energies
I    = ceil((1:(Nx*Nt))'/Nx)                ;
iMix = x.*iL(I) + (1-x).*iG(I)              ;
dMix = 1./(x./delL(I) + (1-x).*1./delG(I))  ;

% %   Back calculate taus
% figure(1);
% h = scatter(dMix,iMix,150,CriticalTemperature()./tau(I),'filled');
% h.Parent.XScale = 'log';
% h.Marker = 'o';
% grid('on');
% box('on');
% axis([1E-5,10,-0.5,9]);


tic;
[~,tauB,~,~,~] = SaturationStateGivenMixedRhoIRRND(dMix,iMix);
toc;


%   Plot
figure(2);
    semilogy((1:numel(tauB))',abs(tau(I) - tauB));

    
figure(3);
h = plot3(dMix,iMix,CriticalTemperature()./tauB,'o');
h.Parent.XScale = 'log';
h.Marker = 'o';
grid('on');
box('on');
axis([1E-5,10,-0.5,9,250,650]);
h.Parent.FontName = 'CMU Concrete';
h.Parent.FontSize = 12;
xl = xlabel('Dimensionless Density [-]','FontName','CMU Concrete','FontSize',14);
xl.Position = [0.0100,-1.1,-1.0000];
yl = ylabel('Dimensionless Internal Energy [-]','FontName','CMU Concrete','FontSize',14);
% yl.Position = [0.0000,4.2500,+1.0000];
c = colorbar('FontName','CMU Concrete','FontSize',12);
c.Label.String = 'Temperature [K]';
c.Label.Rotation = -90;
c.Label.Position = [2.85,460.1277,0];
