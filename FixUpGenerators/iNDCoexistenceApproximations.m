clc;

% clear();
% Tc  = CriticalTemperature();
% tau = Tc./[linspace(273.16,640,1E3),linspace(640.001,Tc,2E3)]';
% [P,delL,delG] = SaturationStateGivenTauRRND(tau);
% iL = InternalEnergyOneRND(delL,tau);
% iG = InternalEnergyOneRND(delG,tau);


% %   Region 1
% mask     = [1,10,300,301,644,645];
% logDelG  = log(delG(mask));
% iGfit    = iG(mask);
% iGfit([2,3,4]) = 1.00008*iGfit([2,3,4]);
% [p,S,mu] = polyfit(logDelG,iGfit,5);
% mask  = mask(1):mask(end);
% delG1 = delG(mask);
% iG1   = polyval(p,log(delG1),S,mu);
% % semilogx(delG(mask),iG(mask),delG(mask),iGfit);
% % legend('Exact','Fit');
% % all(iGfit >= iG(mask))
% 
% %   Region 2
% mask     = [645,646,900,1000,1800,2990,numel(delG)];
% delGfit  = delG(mask);
% iGfit    = iG(mask);
% iGfit([3,4,5]) = 1.0020*iGfit([3,4,5]);
% iGfit(end-1) = 1.0005*iGfit(end-1);
% iGfit(end) = 1.000000000000001*iGfit(end);
% [p,S,mu] = polyfit(delGfit,iGfit,6);
% mask  = mask(1):mask(end);
% delG2 = delG(mask);
% iG2   = polyval(p,delG2,S,mu);
% % plot(delG(mask(1):mask(end)),iG(mask(1):mask(end)),delG2,iG2);
% % legend('Exact','Fit');
% % Show((iG2 - iG(mask))');
% % all(iG2 >= iG(mask))
% 
% semilogx(delG,iG,[delG1;delG2],[iG1;iG2]);
% all([iG1;iG2(2:end)] >= iG)


% %   Region 1
% mask  = [350,675,2000,3000];
% iLfit = iL(mask);
% [p,S,mu] = polyfit(delL(mask),iLfit,3);
% mask = mask(1):mask(end);
% delL1 = delL(mask);
% iL1   = polyval(p,delL1,S,mu);
% shift = 0.02480;
% plot(delL(mask),iL(mask),delL1,iL1+shift);
% Show([mask;(iL1 - iL(mask))']);
% all((iL1+shift) >= iL(mask))


%   Region 2
mask  = [12,50,100,200,350];
iLfit = iL(mask);
[p,S,mu] = polyfit(delL(mask).^2,iLfit,3);
mask = mask(1):mask(end);
delL1 = delL(mask);
iL1   = polyval(p,delL1.^2,S,mu);
shift = 0*0.02;
plot(delL(mask),iL(mask),delL1,iL1+shift);
Show([mask;(iL1 - iL(mask))']);
all((iL1+shift) >= iL(mask))
