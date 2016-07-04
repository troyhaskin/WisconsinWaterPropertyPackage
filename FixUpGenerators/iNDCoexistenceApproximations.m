clc();
clear();



% ============================================= Old fit
%   Pre-stored values
Tc          = CriticalTemperature();
tau         = Tc./[linspace(273.16,280,600),linspace(280.01,640,3E3),linspace(640.001,Tc,2E3)]';
deLmin      = 1.0;
[~,delLmax] = saturableDeltas();
[~,iLAtMax] = iNDsAtSaturableDeltas();
iLmax       = CriticalInternalEnergyND();

%   Calculated values
[P,delL,delG] = SaturationStateGivenTauRRND(tau);
iL = InternalEnergyOneRND(delL,tau);
iG = InternalEnergyOneRND(delG,tau);

mask           = [1,210,450,501,837,838]    ;
logDelG        = log(delG(mask))            ;
iGfit          = iG(mask)                   ;
iGfit([2,3,4]) = 1.0006*iGfit([2,3,4])      ;
iGfit(6)       = 1.000000000000001*iGfit(6) ;
[p,S,mu]       = polyfit(logDelG,iGfit,5)   ;
%   Contract
mask  = mask(1):mask(end);
delG1 = delG(mask);
iG1   = polyval(p,log(delG1),S,mu);

% ============================================= Old fit


%%   New region 1
%
%   Log-Polynomial
order      = 10;
dels       = [delG(1),delG(5500)];
d          = ChebyshevSpace((log(dels(1))),(log(dels(2))),200*order);
[p,~,mu]   = polyfit((-d).^0.1,InternalEnergySaturated(exp(d)),order);
delta      = exp(linspace(log(dels(1)),log(dels(2)),1E3));
iNDexact   = InternalEnergySaturated(delta);
iNDfit     = polyval(p,(-log(delta)).^0.1,[],mu);
while not(all(iNDfit > iNDexact))
    underShoot = min(iNDexact(:) - iNDfit(:));
    p(end)     = p(end) + abs(underShoot);
    iNDfit     = polyval(p,(-log(delta)).^0.1,[],mu);
end
figure(1);
semilogx(delta,iNDexact-iNDfit,'--');
disp(max(iNDfit - iNDexact));
disp(min(iNDfit - iNDexact));
%
disp('Region 1:');
disp('delBound = [');
Show(dels(:).');
disp('];');
disp('p = [');
Show(p);
disp('];');
disp('mu = [');
Show(mu);
disp('];');
%
deltaAll = delta    ;
iNDAll   = iNDfit   ;



%%   New region 2
%
%   Regular Polynomial
disp('Region 2:');
order    = 8;
dels     = [delG(5500),delL(3500)];
d        = ChebyshevSpace((dels(1)),(dels(2)),1E3);
[p,~,mu] = polyfit(log(d),InternalEnergySaturated((d)),order);
delta    = (linspace((dels(1)),(dels(2)),5E4));
iNDexact = InternalEnergySaturated(delta);
iNDfit   = polyval(p,log(delta),[],mu);
while not(all(iNDfit > iNDexact))
    underShoot = min(iNDexact(:) - iNDfit(:));
    p(end)     = p(end) + abs(underShoot);
    iNDfit     = polyval(p,log(delta),[],mu);
end
figure(2);
plot(delta,iNDexact-iNDfit,'--');
disp(max(iNDfit - iNDexact));
disp(min(iNDfit - iNDexact));
%
disp('Region 2:');
disp('delBound = [');
Show(dels(:).');
disp('];');
disp('p = [');
Show(p);
disp('];');
disp('mu = [');
Show(mu);
disp('];');
%
deltaAll = [deltaAll,delta] ;
iNDAll   = [iNDAll,iNDfit];



%%   New region 3 (from near critical to delLmax)
%
%   Regular Polynomial
disp('Region 3:');
order    = 6;
n        = 3500;
dels     = delL(350:n);
iNDexact = iL(350:n);
[p,~,mu] = polyfit(sqrt(delLmax^2- dels.^2),iNDexact,order);
iNDfit   = polyval(p,sqrt(delLmax^2- dels.^2),[],mu);
while not(all(iNDfit > iNDexact))
    underShoot = min(iNDexact - iNDfit);
    p(end)     = p(end) + abs(underShoot);
    iNDfit     = polyval(p,sqrt(delLmax^2- dels.^2),[],mu);
end
figure(3);
plot(dels,iNDexact,dels,iNDfit,'--')
disp(max(iNDfit - iNDexact));
disp(min(iNDfit - iNDexact));
%
disp('Region 3:');
disp('delBound = [');
Show(dels([1,end]).');
disp('];');
disp('p = [');
Show(p);
disp('];');
disp('mu = [');
Show(mu);
disp('];');
%
deltaAll = [deltaAll,dels(end:-1:1).'] ;
iNDAll   = [iNDAll,iNDfit(end:-1:1).'];





%%   New region 4 (from triple to delLmax)
%
%   Regular Polynomial
disp('Region 4:');
order    = 2;
dels     = delL(1:350);
iNDexact = iL(1:350);
[p,~,mu] = polyfit(sqrt(delLmax^2- dels.^2),iNDexact,order);
iNDfit   = polyval(p,sqrt(delLmax^2- dels.^2),[],mu);
while not(all(iNDfit < iNDexact))
    underShoot = max(iNDexact - iNDfit);
    p(end)     = p(end) - abs(underShoot);
    iNDfit     = polyval(p,sqrt(delLmax^2- dels.^2),[],mu);
end
figure(4);
plot(dels,iNDexact,dels,iNDfit,'--')
disp(max(iNDfit - iNDexact));
disp(min(iNDfit - iNDexact));
%
disp('Region 4:');
disp('delBound = [');
Show(dels([1,end]).');
disp('];');
disp('p = [');
Show(p);
disp('];');
disp('mu = [');
Show(mu);
disp('];');
%
deltaAll = [deltaAll,dels(end:-1:1).'] ;
iNDAll   = [iNDAll,iNDfit(end:-1:1).'];


%%

figure(5);
semilogx([delG;delL(end:-1:1)],[iG;iL(end:-1:1)],deltaAll,iNDAll,'--');


% f         = @(delta,c) sum(bsxfun(@times,c(1:end/2).',bsxfun(@power,(log(delta)-mu(1))/mu(2),c((end/2+1):end).')),1);
% objective = @(c) integral(@(delta) abs(f(delta,c) - InternalEnergySaturated(delta)),delG1(1),delG1(end));
% c = fminsearch(objective,c,optimset('Display','iter','MaxIter',30));


