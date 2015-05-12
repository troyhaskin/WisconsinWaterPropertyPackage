clc();

%%
% ===================================================================================== %
%                                        Vapor Dome                                     %
% ===================================================================================== %
clear();
%   Pre-stored values
Tc          = CriticalTemperature();
tau         = Tc./[linspace(273.16,280,200),linspace(280.01,640,1E3),linspace(640.001,Tc,2E3)]';
deLmin      = 1.0;
[~,delLmax] = saturableDeltas();
[~,iLAtMax] = iNDsAtSaturableDeltas();
iLmax       = CriticalInternalEnergyND();

%   Calculated values
[P,delL,delG] = SaturationStateGivenTauRRND(tau);
iL = InternalEnergyOneRND(delL,tau);
iG = InternalEnergyOneRND(delG,tau);



%%
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
fprintf('VaporDome: Gas: Region 1\n');
    fprintf('\tDensity Limits:\n');
        fprintf('\t\t%+23.16E\n',delG(mask([1,end])));
    fprintf('\tCorrelation:\n\t\t%s\n','log-polynomial');
    fprintf('\tp:\n');
        fprintf('\t\t%+23.16E\n',p);
    fprintf('\tmu:\n');
        fprintf('\t\t%+23.16E\n',mu);
    fprintf('\tAll bounded: ')
        if all(iG1 >= iG(mask));
            fprintf('True\n');
        else
            fprintf('False\n');
        end

        
        
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
    fprintf('\tCorrelation:\n\t\t%s\n','log-polynomial');
    fprintf('\tp:\n');
        fprintf('\t\t%+23.16E\n',p);
    fprintf('\tmu:\n');
        fprintf('\t\t%+23.16E\n',mu);
    fprintf('\tAll bounded: ')
        if all(iG2 >= iG(mask));
            fprintf('True\n');
        else
            fprintf('False\n');
        end

        
        
fprintf('\n------------------------------------------------\n');



%%
% =================================================== %
%                    Liquid Region                    %
% =================================================== %

% --------------------------- %
%           Region 1          %
% --------------------------- %
%   Set-up
mask  = 600 : numel(delL);
%   Determine optimal coefficients
[p,S,mu] = polyfit(delL(mask),iL(mask),5);
%   Contract
mask = mask(1):mask(end);
delL1 = delL(mask);
%   SHIFT
shift = +6.1045559161394536E-03;
iL1   = polyval(p,delL1,S,mu) + shift;
%   Print Information
fprintf('VaporDome: Liquid: Region 1\n');
    fprintf('\tDensity Limits:\n');
        fprintf('\t\t%+23.16E\n',delL(mask([1,end])));
    fprintf('\tCorrelation:\n\t\t%s\n','shifted polynomial');
    fprintf('\tp:\n');
        fprintf('\t\t%+23.16E\n',p);
    fprintf('\tmu:\n');
        fprintf('\t\t%+23.16E\n',mu);
    fprintf('\tShoft:\n');
        fprintf('\t\t%+23.16E\n',shift);
    fprintf('\tAll bounded: ')
        if all(iL1 >= iL(mask));
            fprintf('True\n');
        else
            fprintf('False\n');
        end

        
        
fprintf('\n------------------------------------------------\n');



%%
% --------------------------- %
%           Region 2          %
% --------------------------- %
%   Determine fit coefficients
mask  = 300:600;
[p,S,mu] = polyfit(delL(mask),iL(mask),5);
%   Contract
mask = mask(1):mask(end);
delL2 = delL(mask);
%   SHIFT
shift = +7.4979900206617600E-04;
iL2   = polyval(p,delL2,S,mu) + shift;
%   Print Information
fprintf('VaporDome: Liquid: Region 2\n');
    fprintf('\tDensity Limits:\n');
        fprintf('\t\t%+23.16E\n',delL(mask([1,end])));
    fprintf('\tCorrelation:\n\t\t%s\n','shifted polynomial');
    fprintf('\tp:\n');
        fprintf('\t\t%+23.16E\n',p);
    fprintf('\tmu:\n');
        fprintf('\t\t%+23.16E\n',mu);
    fprintf('\tShoft:\n');
        fprintf('\t\t%+23.16E\n',shift);
    fprintf('\tAll bounded: ')
        if all(iL2 >= iL(mask));
            fprintf('True\n');
        else
            fprintf('False\n');
        end

        
        
fprintf('\n------------------------------------------------\n');

        
%%
% --------------------------- %
%           Region 3          %
% --------------------------- %
%   Determine fit coefficients
mask  = [117,200,300];
mask = mask(1):mask(end);
correlation = @(c,x) c(2)*(delLmax - x).^c(1).*exp(c(3)*(delLmax - x)) + iLAtMax;
objective   = @(c) norm(correlation(c,delL(mask)) - iL(mask),Inf);
c = [+5.1770062815506157E-01;+3.2967387240093049E+00;+2.1813809229584362E+00];
c = fminsearch(objective,c);
Show(c);
%   Contract
delL3 = delL(mask);
shift = +5.9565801837926924E-04;
iL3  = correlation(c,delL3) + shift;
%   Print
fprintf('VaporDome: Liquid: Region 3\n');
    fprintf('\tDensity Limits:\n');
        fprintf('\t\t%+23.16E\n',delL(mask([1,end])));
    fprintf('\tCorrelation:\n\t\t%s\n',[func2str(correlation),' + shift']);
    fprintf('\tc:\n');
        fprintf('\t\t%+23.16E\n',c);
    fprintf('\tshift:\n');
        fprintf('\t\t%+23.16E\n',shift);
    fprintf('\tAll bounded: ')
        if all(iL3 >= iL(mask));
            fprintf('True\n');
        else
            fprintf('False\n');
        end

        
        
fprintf('\n------------------------------------------------\n');   


%%
% --------------------------- %
%           Region 4          %
% --------------------------- %
%   Determine fit coefficients
mask  = [1,117];
mask = mask(1):mask(end);
correlation = @(c,x) c(2)*(delLmax - x).^c(1).*exp(c(3)*(delLmax - x))+iLAtMax;
objective   = @(c) norm(correlation(c,delL(mask)) - iL(mask),2);
c = [+4.9763662117192797E-01;-2.7486285142106492E+00;-2.2862334141468203E+01];
c = fminsearch(objective,c);
%   Corrections to ensure over-prediction
shift = -1.1372162114144102E-05;
%   Contract
delL4 = delL(mask);
iL4  = correlation(c,delL4) + shift;
%   Print
fprintf('VaporDome: Liquid: Region 4\n');
    fprintf('\tDensity Limits:\n');
        fprintf('\t\t%+23.16E\n',delL(mask([1,end])));
    fprintf('\tCorrelation:\n\t\t%s\n',[func2str(correlation),' + shift']);
    fprintf('\tc:\n');
        fprintf('\t\t%+23.16E\n',c);
    fprintf('\tshift:\n');
        fprintf('\t\t%+23.16E\n',shift);
    fprintf('\tAll bounded: ')
        if all(iL4 <= iL(mask));
            fprintf('True\n');
        else
            fprintf('False\n');
        end



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
  




