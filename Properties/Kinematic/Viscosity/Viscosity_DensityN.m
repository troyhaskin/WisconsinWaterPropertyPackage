function mu_rhoN = Viscosity_DensityN(rho,T,PhaseCheck,Eps)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    if (nargin < 4) || isempty(PhaseCheck)
        Eps = DefaultDerivativeDelta();
    end
    
    if length(T) > 1
        PropertyHandle = @(rhoVar) Viscosity(rhoVar,[T,T],PhaseCheck);
    else
        PropertyHandle = @(rhoVar) Viscosity(rhoVar,T,PhaseCheck);
    end
    
    mu_rhoN = PointWiseCentralDifference(PropertyHandle,rho,Eps);
    mu_rhoN = RestoreShape(mu_rhoN,GreatestProduct(SizeRho,SizeT));
end