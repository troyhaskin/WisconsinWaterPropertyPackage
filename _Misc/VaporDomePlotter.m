clc;

Font = 'CMU Serif';

Tt  = TriplePointTemperature()  ;
Tc  = CriticalTemperature()     ;
Tsat = linspace(Tt,Tc,1E3);
[Psat,rhol,rhog] = SaturationStateGivenTsat(Tsat);
vl   = 1./rhol;
vg   = 1./rhog;
Psat = Psat/1E6;
Pt  = min(Psat);

ii  = 875;
rho = [linspace(0.005,74,400),linspace(74.0001,648.57999,50),linspace(648.58,997,400)];
v   = 1./rho;
T   = Tsat(ii);
P   = Pressure(rho,T);
P   = P /1E6;

vA = v(479);                            PA = P(479)     ;
vB = vl(ii);                            PB = Psat(ii)   ;
vE = vg(ii);                            PE = PB         ;
vC = (1-0.10)*vB + 0.10*vE;             PC = PB         ;
vD = (1-0.45)*vB + 0.45*vE;             PD = PB         ;
vF = v(200);                            PF = P(200)     ;


semilogx(vl,Psat,'k',vg,Psat,'k','LineWidth',2,'Color',GetColor('grey')); % vapor dome
hold('on');
semilogx(1./rho,P,'k','LineWidth',2); % freezing line
plot(vA,PA,'rs','MarkerFaceColor','r','MarkerSize',8);
plot(vB,PB,'rs','MarkerFaceColor','r','MarkerSize',8);
plot(vC,PC,'rs','MarkerFaceColor','r','MarkerSize',8);
plot(vD,PD,'rs','MarkerFaceColor','r','MarkerSize',8);
plot(vE,PE,'rs','MarkerFaceColor','r','MarkerSize',8);
plot(vF,PF,'rs','MarkerFaceColor','r','MarkerSize',8);
grid('on');
box('on');
xlabel('Volume [m^3/kg]','FontSize',12,'FontName',Font);
ylabel('Pressure [MPa]' ,'FontSize',12,'FontName',Font);
axis([5E-4,1E0,Pt,25]);

text(0.78*vA,0.99*PA,'A','FontSize',12,'FontName',Font);
text(1.01*vB,0.92*PB,'B','FontSize',12,'FontName',Font);
text(1.00*vC,0.92*PC,'C','FontSize',12,'FontName',Font,'HorizontalAlignment','center');
text(1.00*vD,0.92*PD,'D','FontSize',12,'FontName',Font,'HorizontalAlignment','center');
text(1.01*vE,0.92*PE,'E','FontSize',12,'FontName',Font,'HorizontalAlignment','right');
text(1.08*vF,1.04*PF,'F','FontSize',12,'FontName',Font);

set(gca,'FontName',Font);
hold('off');