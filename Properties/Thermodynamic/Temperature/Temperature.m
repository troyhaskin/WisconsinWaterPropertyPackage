function [T,state] = Temperature(rho,i,Tguess,PhaseCheck)
    
    %   Handle inputs and balance sizes
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
    
    %   Nondimensionlize and pass to low-level function
    rhoc        = CriticalDensity()                                 ;
    ic          = DimensioningInternalEnergy()                      ;
    Tc          = CriticalTemperature()                             ;
    delta       = rho / rhoc                                        ;
    iND         = i   / ic                                          ;
    [tau,stateND] = TemperatureRRND(delta,iND,Tc./Tguess,PhaseCheck)  ;
    
    %   Re-dimensionlize variables
    T = Tc./tau;
    state.ND_      = stateND                                ;
    state.rho      = rho                                    ;
    state.i        = i                                      ;
    state.T        = T                                      ;
    state.rhoL     = stateND.delL * rhoc                    ;
    state.rhoG     = stateND.delG * rhoc                    ;
    state.P        = stateND.Pnd  * DimensioningPressure()  ;
    state.x        = stateND.x                              ;
    state.isTwoPhi = stateND.isTwoPhi                       ;
    
    %   Restore shape
    T = RestoreShape(T,GreatestProduct(SizeRho,SizeI));

end