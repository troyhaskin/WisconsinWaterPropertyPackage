function IntEnergy = InternalEnergy(rho,T,PhaseCheck)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if nargin < 3
        PhaseCheck = true;
    end
    
    OnePhiHandle = @(delta,tau,Mask) InternalEnergyOneR(delta,tau);
                                 
    IntEnergy = GenericTwoPhiProperty(rho,T,OnePhiHandle,'quality',PhaseCheck);
    IntEnergy = RestoreShape(IntEnergy,GreatestProduct(SizeRho,SizeT));
    
end
