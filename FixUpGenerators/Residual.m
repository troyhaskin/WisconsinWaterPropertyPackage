

taut = TriplePointTemperatureR();

delL = linspace(3.1049,3.105,500);
delG = linspace(1.5074E-05,1.5078E-05,500);

[delLm,delGm] = meshgrid(delL,delG);
delL = delLm(:);
delG = delGm(:);


PhiR    = @(delta) HelmholtzResidual   (delta,taut);
PhiR_d  = @(delta) HelmholtzResidual_d (delta,taut);
PhiR_dd = @(delta) HelmholtzResidual_dd(delta,taut);

Psig1    = PhiR(delL) - PhiR(delG) + log(delL./delG)    ;
Psig1_dl =  PhiR_d(delL) + 1./delL                      ;
Psig1_dg =-(PhiR_d(delG) + 1./delG)                     ;

PdelL = 1 + delL.*PhiR_d(delL);
PdelG = 1 + delG.*PhiR_d(delG);

PdelL_d = PhiR_d(delL) + delL.*PhiR_dd(delL);
PdelG_d = PhiR_d(delG) + delG.*PhiR_dd(delG);

R1 = delG.*Psig1 - PdelL.*(delL-delG);
R2 = delL.*Psig1 - PdelG.*(delL-delG);

R = R1 + R2;

R = reshape(R,500,500);

H = surf(delLm,delGm,R,'EdgeColor','interp');
% set(get(H,'Parent'),'ZScale','log');
