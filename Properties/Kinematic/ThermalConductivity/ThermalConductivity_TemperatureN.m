function lam_TN = ThermalConductivity_TemperatureN(rho,T,PhaseCheck,Eps)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    if (nargin < 4) || isempty(PhaseCheck)
        Eps = DefaultDerivativeDelta();
    end
    
    PropertyHandle = @(Tvar) ThermalConductivity(rho,Tvar,PhaseCheck);
    
    lam_TN = PointWiseCentralDifference(PropertyHandle,T,Eps);
    lam_TN = RestoreShape(lam_TN,GreatestProduct(SizeRho,SizeT));
end