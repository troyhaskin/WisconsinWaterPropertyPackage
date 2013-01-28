function T_iN = Temperature_InternalEnergyN(rho,Int,Tguess,PhaseCheck,Eps)
    
    [rho,SizeRho,Int,SizeInt] = Columnify(rho,Int);
    
    if (nargin < 3)
        Tguess = [];
    end
    
    if (nargin < 4) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    if (nargin < 5) || isempty(Eps)
        Eps = DefaultDerivativeDelta();
    end
    
    PropertyHandle = @(IntVar) Temperature(rho,IntVar,Tguess,PhaseCheck);
    
    T_iN = PointWiseCentralDifference(PropertyHandle,Int,Eps);
    T_iN = RestoreShape(T_iN,GreatestProduct(SizeRho,SizeInt));

end