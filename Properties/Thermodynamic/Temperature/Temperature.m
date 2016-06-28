function [T,state] = Temperature(rho,i,Tguess,PhaseCheck)
    
    if (nargin < 4)
        PhaseCheck = true;
    end
    
    [rho,SizeRho,i,SizeI] = Columnify(rho,i)        ;
    [rho, i]              = BalanceSizes(rho,i)     ;
    
    if (nargin < 3) || isempty(Tguess)
        Tguess = rho*0 + 600;
    else
        [rho, Tguess] = BalanceSizes(rho,Tguess);
    end
        
% 
% For a given mixture density, if the mixture internal energy
% is greater than the saturated internal energy, it is outside the vapor dome;
% therefore, the trouble of solving the two-phase mixture system can avoided.
%

    delta     = rho / CriticalDensity()                     ;
    iND       = i   / DimensioningInternalEnergy()          ;
    [T,state] = TemperatureRND(delta,iND,Tguess,PhaseCheck) ;


    T = RestoreShape(T,GreatestProduct(SizeRho,SizeI));

end