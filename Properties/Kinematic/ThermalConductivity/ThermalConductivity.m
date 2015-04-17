function lam = ThermalConductivity(rho,T,PhaseCheck)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3)
        PhaseCheck = true;
    end
    
    OnePhiHandle    = @(delta,tau,~) ThermalConductivityOneR(delta,tau);
    
    lam = GenericTwoPhiProperty(rho,T,OnePhiHandle,'void',PhaseCheck);
    lam = RestoreShape(lam,GreatestProduct(SizeRho,SizeT));
    
end