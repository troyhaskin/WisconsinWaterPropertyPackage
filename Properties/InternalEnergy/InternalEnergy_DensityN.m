function i_rhoN = InternalEnergy_DensityN(rho,T,PhaseCheck,Eps)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    if (nargin < 4) || isempty(PhaseCheck)
        Eps = DefaultDerivativeDelta();
    end
    
    if length(T) > 1
        PropertyHandle = @(rhoVar) InternalEnergy(rhoVar,[T,T],PhaseCheck);
    else
        PropertyHandle = @(rhoVar) InternalEnergy(rhoVar,T,PhaseCheck);
    end
    
    i_rhoN = PointWiseCentralDifference(PropertyHandle,rho,Eps);
    i_rhoN = RestoreShape(i_rhoN,GreatestProduct(SizeRho,SizeT));

end