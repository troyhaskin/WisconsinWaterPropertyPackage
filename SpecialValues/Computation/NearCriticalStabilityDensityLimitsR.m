function [delL,delG] = NearCriticalStabilityDensityLimitsR()
    %   Any reduced liquid density below and reduced gas density above these values does 
    %   does not robustly converge to a correct saturation state.  Therefore, 
    %   the densities' and temperature's estimation values will be used with no iteration
    %   attempted.
    delL = 1.37 ;   % 1.4951486610567; 
    delG = 0.65 ;   % 0.54962871618379; 
end