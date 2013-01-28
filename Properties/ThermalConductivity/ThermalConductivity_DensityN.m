function lam_rhoN = ThermalConductivity_DensityN(rho,T,PhaseCheck,Eps)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    if (nargin < 4) || isempty(PhaseCheck)
        Eps = DefaultDerivativeDelta().*rho;
    end
    
    PropertyHandle = @(rhoVar) ThermalConductivity(rhoVar,T,PhaseCheck);
    
    lam_rhoN = PointWiseCentralDifference(PropertyHandle,rho,Eps);
    lam_rhoN = RestoreShape(lam_rhoN,GreatestProduct(SizeRho,SizeT));
end