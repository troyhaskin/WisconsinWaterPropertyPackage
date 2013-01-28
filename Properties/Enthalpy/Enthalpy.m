function h = Enthalpy(rho,T,PhaseCheck)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if nargin < 3
        PhaseCheck = true;
    end

    OnePhiHandle = @(delta,tau,Mask) EnthalpyOneR(delta,tau);
    
    h = GenericTwoPhiProperty(rho,T,OnePhiHandle,'quality',PhaseCheck);
    h = RestoreShape(h,GreatestProduct(SizeRho,SizeT));
end