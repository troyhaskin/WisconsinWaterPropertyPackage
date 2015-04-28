function iNDG = InternalEnergySublimationRND(delta)
%
%   Calculates the dimensionless gas internal energy iNDG along the isochore delta
%   along the sublimation locus.
%
    %   Initialize tau
    tau = delta*0 + TriplePointTau();
    
    %   Solve for tau using pressure equilibrium
    r         = @(tau,mask) PressureOneRND(delta(mask),tau)     - PressureSublimateRND(tau)     ;
    drdtau    = @(tau,mask) PressureOneRND_tau(delta(mask),tau) - PressureSublimateRND_tau(tau) ;
    upupdater = @(tau,mask,r) [r./drdtau(tau,mask),abs(r)]  ;
    updater   = @(tau,mask) upupdater(tau,mask,r(tau,mask)) ;
    tau       = NewtonUpdater(updater,tau,1E-15,30) ;
    
    %   Calculate internal energy
    iNDG = InternalEnergyOneRND(delta,tau);
    
end