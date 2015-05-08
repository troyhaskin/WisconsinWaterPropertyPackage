clc();
clear();

%
%   This script generates the approximations for tau (reduced temperature) near the triple
%   point in the so-called Density Multi-Valued Saturation (DMVS) region where a single
%   saturated liquid density has two corresponding reduced temperatures until the maximum
%   saturation density limit is reached (the maximum reduced density can be accessed via
%   the second output of the saturableDeltas() function).
%

[~,delLMax]         = saturableDeltas();
[tauLo,tauMd,tauHi] = dmvsTaus();
tauBot = linspace(tauLo,tauMd,4000)';
tauTop = linspace(tauMd,tauHi,4000)';
[~,delLBot,~] = SaturationStateGivenTauRRND(tauBot);
[~,delLTop,~] = SaturationStateGivenTauRRND(tauTop);


%%
% --------------------------- %
%         Top Region          %
% --------------------------- %
%   Determine fit coefficients
mask        = [1;2;3;10;50;4000];
correlation = @(c,x) c(2)*(delLMax - x).^c(1).*exp(c(3)*(delLMax - x))+tauTop(mask(1));
objective   = @(c) norm(correlation(c,delLTop(mask)) - tauTop(mask),2);
c = [+4.9928358462384270E-01,+1.6773844460498277E+00,-9.2112364895178622E+00];
c = fminsearch(objective,c);
%   Corrections to ensure over-prediction
c(1) = 0.995*c(1);
c(2) = 0.981*c(2);
%   Contract
tauTopF  = correlation(c,delLTop);
Show(c);
%   Print
figure(1);
subplot(2,1,1);
    plot(delLTop,tauTop,delLTop,tauTopF,'--');
subplot(2,1,2);
    semilogy(delLTop,abs(tauTopF - tauTop));
all(tauTopF >= tauTop)



%%
% --------------------------- %
%        Bottom Region        %
% --------------------------- %
%   Determine fit coefficients
mask        = [1;2;3990;3998;3999;4000];
correlation = @(c,x) c(2)*(delLMax - x).^c(1).*exp(c(3)*(delLMax - x))+tauBot(mask(end));
objective   = @(c) norm(correlation(c,delLBot(mask)) - tauBot(mask),2);
c = [+5.0078747086862485E-01,-1.7094702725985575E+00,+5.9770012891167834E+00];
c = fminsearch(objective,c);
%   Corrections to ensure over-prediction
c(1) = 0.99956*c(1);
c(2) = 1.000*c(2);
%   Contract
tauBotF  = correlation(c,delLBot);
Show(c);
%   Print
figure(2);
subplot(2,1,1);
    plot(delLBot,tauBot,delLBot,tauBotF,'--');
subplot(2,1,2);
    semilogy(delLBot,abs(tauBotF - tauBot));
all(tauBot >= tauBotF)

