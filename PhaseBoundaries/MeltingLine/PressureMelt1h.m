function P = PressureMelt1h(T)
    
    % Determine the temperatures in the correct range
    Mask = (T > 251.165) & (T < 273.16);
    
    % If the mask is not empty, do the work
    if not(isempty(Mask))
        
        % Nondimensionalizers
        Tstar = 273.16   ; %[K]
        Pstar = 611.657  ; %[K]
        
        % Non-dimensional temperature
        Tnd = T(Mask) / Tstar;
        
        % Correlation
        Pnd = 1 + 1.19539337E6 * ( 1 - Tnd.^3.0000E0) + ...
                  8.08183159E4 * ( 1 - Tnd.^2.5750E1) + ...
                  3.33826860E3 * ( 1 - Tnd.^1.0375E2);
        
        % Allocation and selective assignment
        P = T * 0;
        P(Mask) = Pnd * Pstar;

    else
        P = [];
    end
end