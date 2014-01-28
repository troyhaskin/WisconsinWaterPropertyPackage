clc;
clear('all');

Nx = 300   ;
Nt = 50    ;

Tt = TriplePointTemperature()   ;
Tc = CriticalTemperature()      ;

T = linspace(Tt,Tc,Nt)'                     ;
x = GeometricSpace(0,1,1.05,Nx)'            ;
[P,rhoL,rhoG] = SaturationStateGivenTsat(T) ;

rho   = zeros(Nx * Nt,1);
Start =  1 : Nx : ( Nx * Nt);
End   = Nx : Nx : ( Nx * Nt);

for k = 1:Nt
    rho(Start(k):End(k)) = 1./(1./rhoL(k) + x .* (1./rhoG(k) - 1./rhoL(k)));
end

T = bsxfun(@times,ones(Nx,1),T');
T = T(:);

i   = InternalEnergy(rho,T);
iND = i ./ DimensioningInternalEnergy();

Handle = scatter(rho,T,20,iND);
set(get(Handle,'Parent'),'XScale','log');



tauKnown = Tc./T                ;
del = rho / CriticalDensity()   ;

tauBacked = SaturationStateGivenRhoIRRND(del,iND);

    
