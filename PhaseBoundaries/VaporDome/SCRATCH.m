clc;
% clear('all');


iND = 0.5;
del = 3;
tau = 2.3354;

R      = @(t) iND - Helmholtz_t (del,t);
dRdtau = @(t)     - Helmholtz_tt(del,t);


while true
    dtau = R(tau)./dRdtau(tau);
    tau  = tau - dtau;
    
    Show([dtau;tau;R(tau)]);
    if any(abs(dtau) < 100*eps())
        break
    end
end

[Pnd,delL,delG] = SaturationStateGivenTauRRND(tau);

plot(iND,delL,'ro',Helmholtz_t (del,tau),del,'ko');
