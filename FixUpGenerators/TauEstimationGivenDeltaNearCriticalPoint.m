clc;
clear('all');

%   Calculations have shown that EstimateDelLDelGTauFromDel() will generate incorrect
%   (i.e., imaginary or non-physical) guess values for tau if del is "near" 1 since
%   the fit's derivative blows up in that neighorhood.  "Near" means
%       o { deltaL = (1,1.000255) , tau = (1,1.0000000000021) }
%       o { deltaG = (0.99897,1)  , tau = (1.0000000001303,1) }
%
%   Therefore, this file generates the "fix-up" polynomial fit used in 
%   EstimateDelLDelGTauFromDel() for delta values near the critical point. 
%
%   The tau values used will have the same upper bound for both deltaL and deltaG 
%   (for simplicity) and will have the value of 1.0001.  This is several orders of magnitude
%   great than the precision of the neighborhood above but was chosen so that the Estimation 
%   algorithim will be robust (i.e., the critical point is completely removed from the 
%   saturation calculation.

Nfit     = 6E6  ;
OrderFit = 6    ;
tauFit = [linspace(1,1.0000000001,10*Nfit),linspace(1.00000000011,1.00001,Nfit)]';

delLfit = EstimateDelLFromTau(tauFit);
delGfit = EstimateDelGFromTau(tauFit);

[pLfit,sLfit,muLfit] = polyfit(delLfit,tauFit,OrderFit);
delLplot = linspace(1,delLfit(end),25);
tauLplot = polyval(pLfit,delLplot,sLfit,muLfit);

[pGfit,sGfit,muGfit] = polyfit(delGfit,tauFit,OrderFit);
delGplot = linspace(1,delGfit(end),25);
tauGplot = polyval(pGfit,delGplot,sGfit,muGfit);

clf();
figure(1);hold('on');
plot(delLfit,tauFit,'k',delLplot,tauLplot,'ro');
plot(delGfit,tauFit,'k',delGplot,tauGplot,'ro');


Show([tauLplot(1),tauLplot(1)-tauFit(1)]','%+25.16E');
Show([tauLplot(end),tauLplot(end)-tauFit(end)]','%+25.16E');
Show([tauGplot(1),tauGplot(1)-tauFit(1)]','%+25.16E');
Show([tauGplot(end),tauGplot(end)-tauFit(end)]','%+25.16E');

hold('off');

