function i_TN = InternalEnergyR_tauN(delta,tau,PhaseCheck,Eps)
    
    [delta,SizeDelta,tau,SizeTau] = Columnify(delta,tau);
    
    if (nargin < 3) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    if (nargin < 4) || isempty(PhaseCheck)
        Eps = DefaultDerivativeDelta();
    end
    
    if length(rho) > 1
        PropertyHandle = @(tauVar) InternalEnergy([delta,delta],tauVar,PhaseCheck);
    else
        PropertyHandle = @(tauVar) InternalEnergy(delta,tauVar,PhaseCheck);
    end
    
    i_TN = PointWiseCentralDifference(PropertyHandle,tau,Eps);
    i_TN = RestoreShape(i_TN,GreatestProduct(SizeDelta,SizeTau));
end