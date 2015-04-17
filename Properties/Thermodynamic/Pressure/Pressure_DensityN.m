function P_rhoN = Pressure_DensityN(rho,T,PhaseCheck,Eps)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    if (nargin < 4) || isempty(PhaseCheck)
        Eps = DefaultDerivativeDelta();
    end
    
    if length(T) > 1
        PropertyHandle = @(rhoVar) Pressure(rhoVar,[T,T],PhaseCheck);
    else
        PropertyHandle = @(rhoVar) Pressure(rhoVar,T,PhaseCheck);
    end
    
    P_rhoN = PointWiseCentralDifference(PropertyHandle,rho,Eps);
    P_rhoN = RestoreShape(P_rhoN,GreatestProduct(SizeRho,SizeT));

end