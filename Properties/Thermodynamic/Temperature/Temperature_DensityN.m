function T_rhoN = Temperature_DensityN(rho,Int,Tguess,PhaseCheck,Eps)
    
    [rho,SizeRho,Int,SizeInt] = Columnify(rho,Int);
    
    if (nargin < 3) || isempty(Tguess)
        Tguess = [];
    end
    
    if (nargin < 4) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    if (nargin < 5) || isempty(Eps)
        Eps = DefaultDerivativeDelta();
    end
    
    PropertyHandle = @(rhoVar) Temperature(rhoVar,Int,Tguess,PhaseCheck);
    
    T_rhoN = PointWiseCentralDifference(PropertyHandle,rho,Eps);
    T_rhoN = RestoreShape(T_rhoN,GreatestProduct(SizeRho,SizeInt));

end