function mu_TN = Viscosity_TemperatureN(rho,T,PhaseCheck,Eps)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    if (nargin < 4) || isempty(PhaseCheck)
        Eps = DefaultDerivativeDelta();
    end
    
    if length(rho) > 1
        PropertyHandle = @(Tvar) Viscosity([rho,rho],Tvar,PhaseCheck);
    else
        PropertyHandle = @(Tvar) Viscosity(rho,Tvar,PhaseCheck);
    end
    
    mu_TN = PointWiseCentralDifference(PropertyHandle,T,Eps);
    mu_TN = RestoreShape(mu_TN,GreatestProduct(SizeRho,SizeT));
end