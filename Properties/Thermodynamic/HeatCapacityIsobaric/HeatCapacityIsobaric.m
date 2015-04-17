function cp = HeatCapacityIsobaric(rho,T,PhaseCheck,TwoPhaseHandle)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3)
        PhaseCheck = true;
    end
    
    if (nargin < 4)
        TwoPhaseHandle = 'Homogeneous Void Fraction';
    end
    
    
    OnePhiHandle    = @(delta,tau,Mask) HeatCapacityIsobaricSingleR(delta,tau);
    
    cp = GenericTwoPhiProperty(rho,T,OnePhiHandle,TwoPhaseHandle,PhaseCheck);
    cp = RestoreShape(cp,GreatestProduct(SizeRho,SizeT));
    
end