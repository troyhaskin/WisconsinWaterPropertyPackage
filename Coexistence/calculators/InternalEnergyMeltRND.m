function iNDL = InternalEnergyMeltRND(delta)
%
%   Calculates the dimensionless liquid internal energy iNDL along the isochore delta
%   along the melt locus for Ice-Ih.
%
    %   Initialize tau
    tau = delta*0 + TriplePointTau();
    
    %   Solve for tau using pressure equilibrium
    r         = @(tau,mask) PressureOneRND(delta(mask),tau)     - PressureMeltIhRND(tau)     ;
    drdtau    = @(tau,mask) PressureOneRND_tau(delta(mask),tau) - PressureMeltIhRND_tau(tau) ;
    upupdater = @(tau,mask,r) [r./drdtau(tau,mask),abs(r)]  ;
    updater   = @(tau,mask) upupdater(tau,mask,r(tau,mask)) ;
    tau       = NewtonUpdater(updater,tau,1E-15,30) ;
    
    %   Calculate internal energy
    iNDL = InternalEnergyOneRND(delta,tau);
    
end