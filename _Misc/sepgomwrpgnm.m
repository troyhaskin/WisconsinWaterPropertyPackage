clc;
% clear('all');
% 
% % % tau  = CriticalTemperature ./ ChebyshevSpace(647,CriticalTemperature,10E3);
TsingleValue = 281.3143;
% T  = ChebyshevSpace(TriplePointTemperature,350,50E3);
% T  = ChebyshevSpace(350,500,100E3);
% T  = ChebyshevSpace(500,625,100E3);
% T  = ChebyshevSpace(625,640,5E4);
T  = ChebyshevSpace(640,CriticalTemperature,1E6);
[Psat,rhol,rhog] = SaturationStateGivenTsat(T);
tau = CriticalTemperature() ./ T;
delG = rhog / CriticalDensity();
delL = rhol / CriticalDensity();


% del = delG;
% [p,~,mu] = polyfit(del,tau,13);
% figure(1);
% plot(del,polyval(p,del,[],mu),'ro',del,tau);
% Show(norm(polyval(p,del,[],mu)-tau,1));
% Show(norm(polyval(p,del,[],mu)-tau,1)/norm(tau,1));
% 
% 
% del = delL;
% [p,~,mu] = polyfit(del,tau,13);
% figure(2);
% plot(del,polyval(p,del,[],mu),'ro',del,tau);
% Show(norm(polyval(p,del,[],mu)-tau,1));
% Show(norm(polyval(p,del,[],mu)-tau,1)/norm(tau,1));

%%
[p,~,mu] = polyfit(delG,delL,6);
figure(3);
plot(delG,polyval(p,delG,[],mu),'ro',delG,delL);
disp(' ');disp('delL vs. delG');
Show(polyval(p,1,[],mu));
Show(norm(polyval(p,delG,[],mu)-delL,1));
Show(norm(polyval(p,delG,[],mu)-delL,1)/norm(delL,1));

%%
[p,~,mu] = polyfit(delL,delG,6);
figure(4);
plot(delL,polyval(p,delL,[],mu),'ro',delL,delG);
disp(' ');disp('delG vs. delL');
Show(polyval(p,1,[],mu));
Show(norm(polyval(p,delL,[],mu)-delG,1));
Show(norm(polyval(p,delL,[],mu)-delG,1)/norm(delG,1));

% Show(p,'%+20.10E')
% disp(' ');
% Show(mu','%+20.10E')

