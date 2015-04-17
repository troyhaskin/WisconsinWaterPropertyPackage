function i_TN = InternalEnergy_TemperatureN(rho,T,PhaseCheck,Eps)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    if (nargin < 4) || isempty(PhaseCheck)
        Eps = DefaultDerivativeDelta();
    end
    
    if length(rho) > 1
        PropertyHandle = @(Tvar) InternalEnergy([rho,rho],Tvar,PhaseCheck);
    else
        PropertyHandle = @(Tvar) InternalEnergy(rho,Tvar,PhaseCheck);
    end
    
    i_TN = PointWiseCentralDifference(PropertyHandle,T,Eps);
    i_TN = RestoreShape(i_TN,GreatestProduct(SizeRho,SizeT));
end