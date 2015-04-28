clc;

%%
% ===================================================================================== %
%                                        Vapor Dome                                     %
% ===================================================================================== %
clear();
Tc  = CriticalTemperature();
tau = Tc./[linspace(273.16,280,200),linspace(280.01,640,1E3),linspace(640.001,Tc,2E3)]';
[P,delL,delG] = SaturationStateGivenTauRRND(tau);
iL = InternalEnergyOneRND(delL,tau);
iG = InternalEnergyOneRND(delG,tau);


% =================================================== %
%                      Gas Region                     %
% =================================================== %
%%
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
fprintf('VaporDome: Gas: Region 1\n');
    fprintf('\tDensity Limits:\n');
        fprintf('\t\t%+23.16E\n',delG(mask([1,end])));
    fprintf('\tCorrelation:\n\t\t%s\n','polynomial');
    fprintf('\tp:\n');
        fprintf('\t\t%+23.16E\n',p);
    fprintf('\tmu:\n');
        fprintf('\t\t%+23.16E\n',mu);

        
        
fprintf('\n------------------------------------------------\n');


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
fprintf('VaporDome: Gas: Region 2\n');
    fprintf('\tDensity Limits:\n');
        fprintf('\t\t%+23.16E\n',delG(mask([1,end])));
    fprintf('\tCorrelation:\n\t\t%s\n','polynomial');
    fprintf('\tp:\n');
        fprintf('\t\t%+23.16E\n',p);
    fprintf('\tmu:\n');
        fprintf('\t\t%+23.16E\n',mu);

        
        
fprintf('\n------------------------------------------------\n');




% =================================================== %
%                    Liquid Region                    %
% =================================================== %
%%
% --------------------------- %
%           Region 1          %
% --------------------------- %
%   Determine fit coefficients
mask  = 1095:numel(delL);
correlation = @(c,x) c(1)*(x - delL(mask(end))).^c(2).*exp(c(3)*( c(4)*x - delL(mask(end)))) + iL(mask(end));
objective   = @(c) norm(correlation(c,delL(mask)) - 1.0000*iL(mask),1);
c = [+4.64,+0.639,-0.04878,+1.513];
c = fminsearch(objective,c);
%   Contract
mask = mask(1):mask(end);
delL1 = delL(mask);
%   SHIFT
shift = 5.53E-3;
iL1   = correlation(c,delL1) + shift;
% plot(delL(mask),iL(mask),delL1,iL1,'--');
%   Print Information
fprintf('VaporDome: Liquid: Region 1\n');
    fprintf('\tDensity Limits:\n');
        fprintf('\t\t%+23.16E\n',delL(mask([1,end])));
    fprintf('\tCorrelation:\n\t\t%s\n',[func2str(correlation),' + shift']);
    fprintf('\tdelL(mask(end)):\n\t\t%+23.16E\n',delL(mask(end)));
    fprintf('\tiL(mask(end)):\n\t\t%+23.16E\n'  ,iL(mask(end)));
    fprintf('\tshift:\n\t\t%+23.16E\n'  ,shift);
    fprintf('\tc:\n');
        fprintf('\t\t%+23.16E\n',c);


fprintf('\n------------------------------------------------\n');


% --------------------------- %
%           Region 2          %
% --------------------------- %
%   Determine fit coefficients
mask  = 204:1095;
correlation = @(c,x) c(1)*(delL(mask(1)) - x).^c(2).*exp(c(3)*(delL(mask(1)) - c(4)*x)) + iL(mask(1));
objective   = @(c) norm(correlation(c,delL(mask)) - iL(mask),2);
c = [+4.6349302859479611E+00,+6.3859197437871751E-01,-4.8779370095661273E-02,+1.5127244173436336E+00];
c = fminsearch(objective,c);
%   Contract
mask = mask(1):mask(end);
delL2 = delL(mask);
%   SHIFT
shift = 1.386E-2;
iL2   = correlation(c,delL2) + shift;
%   Print Information
fprintf('VaporDome: Liquid: Region 2\n');
    fprintf('\tDensity Limits:\n');
        fprintf('\t\t%+23.16E\n',delL(mask([1,end])));
    fprintf('\tCorrelation:\n\t\t%s\n',[func2str(correlation),' + shift']);
    fprintf('\tdelL(mask(end)):\n\t\t%+23.16E\n',delL(mask(1)));
    fprintf('\tiL(mask(end)):\n\t\t%+23.16E\n'  ,iL(mask(1)));
    fprintf('\tshift:\n\t\t%+23.16E\n'  ,shift);
    fprintf('\tc:\n');
        fprintf('\t\t%+23.16E\n',c);

        
        
fprintf('\n------------------------------------------------\n');

        

% --------------------------- %
%           Region 3          %
% --------------------------- %
%   Determine fit coefficients
mask  = [117,200,204];
mask = mask(1):mask(end);
correlation = @(c,x) c(2)*(max(delL(mask)) - x).^c(1).*exp(c(3)*(max(delL) - x))+iL(mask(1));
objective   = @(c) norm(correlation(c,delL(mask)) - iL(mask),2);
c = [0.49,2.67,15];
c = fminsearch(objective,c);
%   Corrections to ensure over-prediction
c(1) = 0.980*c(1);
c(2) = 0.927*c(2);
%   Contract
delL3 = delL(mask);
iL3  = correlation(c,delL3);
%   Print
fprintf('VaporDome: Liquid: Region 3\n');
    fprintf('\tDensity Limits:\n');
        fprintf('\t\t%+23.16E\n',delL(mask([1,end])));
    fprintf('\tCorrelation:\n\t\t%s\n',func2str(correlation));
    fprintf('\t(max(delL(mask)):\n\t\t%+23.16E\n',max(delL(mask)));
    fprintf('\tiL(mask(1)):\n\t\t%+23.16E\n'  ,iL(mask(1)));
    fprintf('\tc:\n');
        fprintf('\t\t%+23.16E\n',c);



fprintf('\n------------------------------------------------\n');      


% --------------------------- %
%           Region 4          %
% --------------------------- %
%   Determine fit coefficients
mask  = [1,117];
mask = mask(1):mask(end);
correlation = @(c,x) c(2)*(max(delL(mask)) - x).^c(1).*exp(c(3)*(max(delL) - x))+iL(mask(end));
objective   = @(c) norm(correlation(c,delL(mask)) - iL(mask),2);
c = [0.49,2.67,15];
c = fminsearch(objective,c);
%   Corrections to ensure over-prediction
c(1) = 0.9998*c(1);
c(2) = 1.000*c(2);
%   Contract
delL4 = delL(mask);
iL4  = correlation(c,delL4);
%   Print
fprintf('VaporDome: Liquid: Region 4\n');
    fprintf('\tDensity Limits:\n');
        fprintf('\t\t%+23.16E\n',delL(mask([1,end])));
    fprintf('\tCorrelation:\n\t\t%s\n',func2str(correlation));
    fprintf('\t(max(delL(mask)):\n\t\t%+23.16E\n',max(delL(mask)));
    fprintf('\tiL(mask(end)):\n\t\t%+23.16E\n'  ,iL(mask(end)));
    fprintf('\tc:\n');
        fprintf('\t\t%+23.16E\n',c);



fprintf('\n------------------------------------------------\n');



% ===================================================================================== %
%                                     Melting Lines                                     %
% ===================================================================================== %
%%
% --------------------------------------------------- %
%                Liquid Melting Line                  %
% --------------------------------------------------- %

%   Melt curve: liquid side
[Pnd,tauM] = MeltingLineRND([5E3,0,0,0,0]);
delL0     = 1E3/CriticalDensity();
upupdater = @(delL,mask,r) [r./PressureOneRND_delta(delL,tauM(mask)),abs(r)];
updater   = @(delL,mask) upupdater(delL,mask,PressureOneRND(delL,tauM(mask)) - Pnd(mask));
delM  = NewtonUpdater(updater,Pnd*0+delL0,1E-15,30);
iNDM  = InternalEnergyOneRND(delM,tauM);

mask  = [1,900,3500,numel(delM)];
iNDMfit = iNDM(mask);
iNDMfit(2) = 0.99*iNDMfit(2);
iNDMfit(3) = 0.996*iNDMfit(3);
[p,S,mu] = polyfit(delM(mask),iNDMfit,3);
delM1 = delM;
iNDM1 = polyval(p,delM1,S,mu);
% plot(delM,iNDM,delM1,iNDM1,delM1(mask),iNDM1(mask),'o');
% Show((iNDM1(1:1) - iNDM(1:1))');
% all(iNDM(1:end) <= iNDM1(1:end));
fprintf('Melt Curve: Liquid Line\n');
    fprintf('\tDensity Limits:\n');
        fprintf('\t\t%+23.16E\n',[min(delM),max(delM)]);
    fprintf('\tCorrelation:\n\t\t%s\n','polynomial');
    fprintf('\tp:\n');
        fprintf('\t\t%+23.16E\n',p);
    fprintf('\tmu:\n');
        fprintf('\t\t%+23.16E\n',mu);

        
        
fprintf('\n------------------------------------------------\n');



% --------------------------------------------------- %
%                   Gas Melting Line                  %
% --------------------------------------------------- %

%   Melt curve: liquid side
[Pnd,tau] = SublimationLineRND(10E3);
mask      = CriticalTemperature()./tau >= 200;
tau       = tau(mask);
Pnd       = Pnd(mask);
del0      = Pnd .* tau;
upupdater = @(del,mask,r) [r./PressureOneRND_delta(del,tau(mask)),abs(r)];
updater   = @(del,mask) upupdater(del,mask,PressureOneRND(del,tau(mask)) - Pnd(mask));
delS       = NewtonUpdater(updater,del0,1E-15,30);
iNDS       = InternalEnergyOneRND(delS,tau);

mask      = [1,99,1200,2600,numel(delS)] ;
iNDfit    = iNDS(mask)               ;
iNDfit([1,5]) = 1.000000000000001*iNDfit([1,5]);
iNDfit(2) = 1.00002*iNDfit(2)          ;
iNDfit(3) = 1.00005*iNDfit(3)         ;
iNDfit(4) = 1.00001*iNDfit(4)         ;
[p,S,mu]  = polyfit(log(delS(mask)),iNDfit,4);
delS1       = delS;
iNDS1       = polyval(p,log(delS1),S,mu);
% semilogx(delS,iNDS,delS1,iNDS1,delS(mask),iNDS(mask),'o');
% Show((iNDS1(1:end) - iNDS(1:end))');
% all(iNDS(1:end) <= iNDS1(1:end));
fprintf('Melt Curve: Gas Line\n');
    fprintf('\tDensity Limits:\n');
        fprintf('\t\t%+23.16E\n',[min(delS),max(delS)]);
    fprintf('\tCorrelation:\n\t\t%s\n','polynomial');
    fprintf('\tp:\n');
        fprintf('\t\t%+23.16E\n',p);
    fprintf('\tmu:\n');
        fprintf('\t\t%+23.16E\n',mu);

        
        
fprintf('\n------------------------------------------------\n');




% =================================================== %
%                         Plot                        % 
% =================================================== %
%%
figure(1);

subplot(2,1,1);
    semilogx(delG,iG,[delG1;delG2(2:end)],[iG1;iG2(2:end)]);
    xlabel('\delta_g')
    ylabel('i_{nd}');
    legend('Exact','Fit','Location','southwest');
subplot(2,1,2);
    loglog(delG,[iG1;iG2(2:end)] - iG);
    xlabel('\delta_g')
    ylabel('\Delta{i_{nd}}');

figure(2);
subplot(2,1,1);
    plot(delL,iL,[delL4;delL3(2:end);delL2(2:end);delL1(2:end)],[iL4;iL3(2:end);iL2(2:end);iL1(2:end)]);
    xlabel('\delta_l')
    ylabel('i_{nd}');
    legend('Exact','Fit','Location','southwest');
subplot(2,1,2);
    semilogy([delL3;delL2(2:end);delL1(2:end)],[iL3;iL2(2:end);iL1(2:end)] - iL(117:end));
    xlabel('\delta_l')
    ylabel('\Delta{i_{nd}}');
    
figure(3);
subplot(2,1,1);
    plot(delM,iNDM,delM1,iNDM1);
    xlabel('\delta_l')
    ylabel('i_{nd}');
    legend('Exact','Fit','Location','southwest');
subplot(2,1,2);
    semilogy(delM,iNDM1 - iNDM);
    xlabel('\delta_l')
    ylabel('\Delta{i_{nd}}');
    
    
figure(4);
subplot(2,1,1);
    plot(delS,iNDS,delS1,iNDS1);
    xlabel('\delta_l')
    ylabel('i_{nd}');
    legend('Exact','Fit','Location','southwest');
subplot(2,1,2);
    semilogy(delS,iNDS1 - iNDS);
    xlabel('\delta_l')
    ylabel('\Delta{i_{nd}}');
  




