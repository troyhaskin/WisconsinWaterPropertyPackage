function P_TN = Pressure_TemperatureN(rho,T,PhaseCheck,Eps)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    if (nargin < 4) || isempty(PhaseCheck)
        Eps = DefaultDerivativeDelta();
    end
    
    PropertyHandle = @(Tvar) Pressure(rho,Tvar,PhaseCheck);
    
    P_TN = PointWiseCentralDifference(PropertyHandle,T,Eps);
    P_TN = RestoreShape(P_TN,GreatestProduct(SizeRho,SizeT));
end