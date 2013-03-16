% clc;
% clear('all');

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

% Fitting parameters
Nfit     = 6.5E5  ;
OrderFit = 6;

% taus used to sample the guessing correlation
tauSamp = ChebyshevSpace(1,1.00001,Nfit);

% liquid and gas deltas sampled from the guessing correlation
delLsamp = EstimateDelLFromTau(tauSamp);
delGsamp = EstimateDelGFromTau(tauSamp);


% polynomial fit for tau as a function of liquid delta
[pLfit,sLfit,muLfit] = polyfit(delLsamp,tauSamp,OrderFit);
delLFit = linspace(1,delLsamp(end),25);
tauLFit = polyval(pLfit,delLFit,sLfit,muLfit);

% polynomial fit for tau as a function of gas delta
[pGfit,sGfit,muGfit] = polyfit(delGsamp,tauSamp,OrderFit);
delGFit = linspace(1,delGsamp(end),25);
tauGFit = polyval(pGfit,delGFit,sGfit,muGfit);

% Sample values for Fit points
tauSampLCompare = interp1(delLsamp,tauSamp,delLFit,'pchip');
tauSampGCompare = interp1(delGsamp,tauSamp,delGFit,'pchip');


% Plot correlation and fit values
figure(1);
plot(delLsamp,tauSamp,'k',delLFit,tauLFit,'ro',...
     delGsamp,tauSamp,'k',delGFit,tauGFit,'ro');
box('on');
grid('on');
xlabel('delta [-]');
ylabel('tau [-]');

figure(2);
semilogy(delLFit,abs(tauSampLCompare - tauLFit),...
         delGFit,abs(tauSampGCompare - tauGFit));
box('on');
grid('on');
xlabel('delta [-]');
ylabel('Delta tau [-]');



% Norms and Statistics
fprintf('                Summary of Fit Statistics              \n');
fprintf('=======================================================\n');
fprintf('\nOrder of polynomial: %G\n\n',OrderFit);
% Reduced Liquid Fit Statistics
fprintf('Reduced Liquid Fit:\n');
fprintf('    Critical tau:\n');
fprintf('        Value:          %+23.16E\n',tauLFit(1));
fprintf('        Absolute Error: %+23.16E\n\n',tauLFit(1) - 1);
fprintf('    Compared to Correlation:\n');
fprintf('        L_1   Norm:      %+23.16E\n',norm(tauSampLCompare - tauLFit,1));
fprintf('        L_Inf Norm:      %+23.16E\n\n',norm(tauSampLCompare - tauLFit,Inf));
% Reduced Liquid Fit Statistics
fprintf('Reduced Gas Fit:\n');
fprintf('    Critical tau:\n');
fprintf('        Value:          %+23.16E\n',tauGFit(1));
fprintf('        Absolute Error: %+23.16E\n\n',tauGFit(1) - 1);
fprintf('    Compared to Correlation:\n');
fprintf('        L_1   Norm:     %+23.16E\n',norm(tauSampGCompare - tauGFit,1));
fprintf('        L_Inf Norm:     %+23.16E\n',norm(tauSampGCompare - tauGFit,Inf));


