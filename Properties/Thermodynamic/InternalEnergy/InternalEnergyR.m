function IntEnergy = InternalEnergyR(delta,tau,PhaseCheck)
    
    [delta,SizeDelta,tau,SizeTau] = Columnify(delta,tau);
    
    if nargin < 3
        PhaseCheck = true;
    end
    
    OnePhiHandle = @(delta,tau,Mask) InternalEnergyOneR(delta,tau);
                                 
    IntEnergy = GenericTwoPhiPropertyR(delta,tau,OnePhiHandle,'quality',PhaseCheck);
    IntEnergy = RestoreShape(IntEnergy,GreatestProduct(SizeDelta,SizeTau));
    
end
