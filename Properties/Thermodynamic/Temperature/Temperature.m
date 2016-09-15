function [T,varargout] = Temperature(rho,i,Tguess,PhaseCheck)
    
    %   Handle inputs and balance sizes
    if (nargin < 4)
        PhaseCheck = true;
    end
    
    [rho,SizeRho,i,SizeI] = Columnify(rho,i)        ;
    [rho, i]              = BalanceSizes(rho,i)     ;
    
    %   Load constants
    rhoc = CriticalDensity()            ;
    ic   = DimensioningInternalEnergy() ;
    Tc   = CriticalTemperature()        ;
    
    %   Get good guess
    if (nargin < 3) || isempty(Tguess)
        Tguess = (rho < rhoc)*700 + (rho >= rhoc)*500;
    else
        [rho, Tguess] = BalanceSizes(rho,Tguess);
    end
    mask = isnan(Tguess);
    if any(mask)
        Tguess(mask) = (rho(mask) < rhoc)*700 + (rho(mask) >= rhoc)*500;
    end
    
   
    %   Nondimensionlize and pass to low-level function
    delta       = rho / rhoc                                        ;
    iND         = i   / ic                                          ;
    [tau,stateND] = TemperatureRRND(delta,iND,Tc./Tguess,PhaseCheck)  ;
    
    %   Re-dimensionlize variable and restore shape
    T = Tc./real(tau)                                   ;
    T = RestoreShape(T,GreatestProduct(SizeRho,SizeI))  ;
    
    if (nargout > 1)
        restore        = @(x) RestoreShape(x,GreatestProduct(SizeRho,SizeI))                ;
        state.ND_      = stateND                                                            ;
        state.rho      = restore(rho)                                                       ;
        state.i        = restore(i)                                                         ;
        state.T        = T                                                                  ;
        state.P        = restore(stateND.Pnd  * DimensioningPressure())                     ;
        state.isTwoPhi = restore(stateND.isTwoPhi)                                          ;
        state.x        = restore(stateND.x.*stateND.isTwoPhi - not(stateND.isTwoPhi)*100)   ;
        state.rhoL     = restore(stateND.delL) * rhoc                                       ;
        state.rhoG     = restore(stateND.delG) * rhoc                                       ;
        varargout{1}   = state                                                              ;
    end

end