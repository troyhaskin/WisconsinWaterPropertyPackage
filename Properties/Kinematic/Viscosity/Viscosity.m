function mu = Viscosity(rho,T,PhaseCheck,state)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    OnePhiHandle    = @(delta,tau,Mask) ViscosityOneR(delta,tau);
    
    mu = GenericTwoPhiProperty(rho,T,OnePhiHandle,'void',PhaseCheck,state);
    mu = RestoreShape(mu,GreatestProduct(SizeRho,SizeT));
    
end