clc;

% clear();
% Tc  = CriticalTemperature();
% tau = Tc./[linspace(273.16,280,200),linspace(280.01,640,1E3),linspace(640.001,Tc,2E3)]';
% [P,delL,delG] = SaturationStateGivenTauRRND(tau);
% iL = InternalEnergyOneRND(delL,tau);
% iG = InternalEnergyOneRND(delG,tau);


% =================================================== %
%                      Gas Region                     %
% =================================================== %

% --------------------------- %
%           Region 1          %
% --------------------------- %
%   Indices of interest
mask     = [1,210,450,501,837,838];
logDelG  = log(delG(mask));
iGfit    = iG(mask);
%   Correction to ensure over-prediction
iGfit([2,3,4]) = 1.0006*iGfit([2,3,4]);
iGfit(6) = 1.000000000000001*iGfit(6);
%   Determine coefficients
[p,S,mu] = polyfit(logDelG,iGfit,5);
%   Contract
mask  = mask(1):mask(end);
delG1 = delG(mask);
iG1   = polyval(p,log(delG1),S,mu);


% --------------------------- %
%           Region 2          %
% --------------------------- %
%   Indices of interest
mask     = [838,839,1093,1193,1993,3183,numel(delG)];
delGfit  = delG(mask);
iGfit    = iG(mask);
%   Correction to ensure over-prediction
iGfit([3,4,5]) = 1.0020*iGfit([3,4,5]);
iGfit(end-1) = 1.0006*iGfit(end-1);
iGfit(end) = 1.000000000000001*iGfit(end);
%   Determine coefficients
[p,S,mu] = polyfit(delGfit,iGfit,6);
%   Contract
mask  = mask(1):mask(end);
delG2 = delG(mask);
iG2   = polyval(p,delG2,S,mu);






% =================================================== %
%                    Liquid Region                    %
% =================================================== %

% --------------------------- %
%           Region 1          %
% --------------------------- %
%   Indices of interest
mask  = [206,207,1000,1800,2800,numel(delL)];
iLfit = iL(mask);
%   Correction to ensure over-prediction
iLfit(1) = 1.000000001*iLfit(1);
iLfit(3) = 1.002*iLfit(3);
iLfit(4) = 1.002*iLfit(4);
iLfit(5) = 1.002*iLfit(5);
iLfit(end) = 1.000000000000001*iLfit(end);
%   Determine coefficients
[p,S,mu] = polyfit(delL(mask),iLfit.^2,5);
%   Contract
mask = mask(1):mask(end);
delL1 = delL(mask);
iL1   = sqrt(polyval(p,delL1,S,mu));


% --------------------------- %
%           Region 2          %
% --------------------------- %
%   Determine fit coefficients
mask  = [117,200,206];
mask = mask(1):mask(end);
correlation = @(c,x) c(2)*(max(delL(mask)) - x).^c(1).*exp(c(3)*(max(delL) - x))+iL(117);
objective   = @(c) norm(correlation(c,delL(mask)) - iL(mask),2);
c = [0.49,2.67,15];
c = fminsearch(objective,c);
%   Corrections to ensure over-prediction
c(1) = 0.980*c(1);
c(2) = 0.927*c(2);
%   Contract
delL2 = delL(mask);
iL2  = correlation(c,delL2);



% --------------------------- %
%           Region 3          %
% --------------------------- %
%   Determine fit coefficients
mask  = [1,117];
mask = mask(1):mask(end);
correlation = @(c,x) c(2)*(max(delL(mask)) - x).^c(1).*exp(c(3)*(max(delL) - x))+iL(117);
objective   = @(c) norm(correlation(c,delL(mask)) - iL(mask),2);
c = [0.49,2.67,15];
c = fminsearch(objective,c);
%   Corrections to ensure over-prediction
c(1) = 0.9998*c(1);
c(2) = 1.000*c(2);
%   Contract
delL3 = delL(mask);
iL3  = correlation(c,delL3);

plot(delL(1:206),iL(1:206),[delL3;delL2],[iL3;iL2]);






% =================================================== %
%               Ice Coexistence Regions               %
% =================================================== %

% %   Melt curve: liquid side
% [Pnd,tauM] = MeltingLineRND([5E3,0,0,0,0]);
% delL0     = 1E3/CriticalDensity();
% upupdater = @(delL,mask,r) [r./PressureOneRND_delta(delL,tauM(mask)),abs(r)];
% updater   = @(delL,mask) upupdater(delL,mask,PressureOneRND(delL,tauM(mask)) - Pnd(mask));
% delM      = NewtonUpdater(updater,Pnd*0+delL0,1E-15,30);
% iNDM      = InternalEnergyOneRND(delM,tauM);
% plot(delL2,iL2,delM,iNDM);


% 
% % =================================================== %
% %                         Plot                        % 
% % =================================================== %
% 
% figure(1);
% 
% subplot(2,2,1);
%     semilogx(delG,iG,[delG1;delG2(2:end)],[iG1;iG2(2:end)]);
%     xlabel('\delta_g')
%     ylabel('i_{nd}');
%     legend('Exact','Legend','Location','southwest');
% subplot(2,2,2);
%     loglog(delG,[iG1;iG2(2:end)] - iG);
%     xlabel('\delta_g')
%     ylabel('\Delta{i_{nd}}');
% 
% 
% subplot(2,2,3);
%     plot(delL,iL,[delL2;delL1(2:end)],[iL2;iL1(2:end)]);
%     xlabel('\delta_l')
%     ylabel('i_{nd}');
%     legend('Exact','Legend','Location','southwest');
% subplot(2,2,4);
%     semilogy([delL2;delL1(2:end)],[iL2;iL1(2:end)] - iL(117:end));
%     xlabel('\delta_l')
%     ylabel('\Delta{i_{nd}}');
% 




