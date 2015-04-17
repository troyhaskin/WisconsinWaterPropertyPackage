function cv = HeatCapacityIsochoric(rho,T,PhaseCheck,TwoPhaseHandle)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    if (nargin < 4)
        TwoPhaseHandle = 'Homogeneous Void Fraction';
    end
    
    OnePhiHandle    = @(delta,tau,~) HeatCapacityIsochoricOneR(delta,tau);
    
    cv = GenericTwoPhiProperty(rho,T,OnePhiHandle,TwoPhaseHandle,PhaseCheck);
    cv = RestoreShape(cv,GreatestProduct(SizeRho,SizeT));
end