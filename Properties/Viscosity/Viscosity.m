function mu = Viscosity(rho,T,PhaseCheck)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3)
        PhaseCheck = true;
    end
    
    OnePhiHandle    = @(delta,tau,Mask) ViscosityOneR(delta,tau);
    
    mu = GenericTwoPhiProperty(rho,T,OnePhiHandle,'void',PhaseCheck);
    mu = RestoreShape(mu,GreatestProduct(SizeRho,SizeT));
    
end