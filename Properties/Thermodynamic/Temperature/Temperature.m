function T = Temperature(rho,i,Tguess,PhaseCheck)
    
    if (nargin < 4)
        PhaseCheck = true;
    end
    
    [rho,SizeRho,i,SizeI] = Columnify(rho,i)        ;
    [rho, i]              = BalanceSizes(rho,i)     ;
    
    if (nargin < 3)
        Tguess = rho*0 + 300                    ;
    else
        [rho, Tguess] = BalanceSizes(rho,Tguess);
    end
        
% 
% For a given mixture density, if the mixture internal energy
% is greater than the saturated internal energy, it is outside the vapor dome;
% therefore, the trouble of solving the two-phase mixture system can avoided.
%

    delta = rho / CriticalDensity()                  ;
    iND   = i   / DimensioningInternalEnergy()       ;
    T     = TemperatureRND(delta,iND,Tguess,PhaseCheck) ;


    T = RestoreShape(T,GreatestProduct(SizeRho,SizeI));

end