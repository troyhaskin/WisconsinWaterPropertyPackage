function P = PressureMeltVI(T)
    
    % Determine the temperatures in the correct range
    Mask = (T > 273.31) & (T < 355.0);
    
    % If the mask is not empty, do the work
    if not(isempty(Mask))
        
        % Nondimensionalizers
        Tstar = 273.31     ; %[K]
        Pstar = 6.3240E+08  ; %[K]
        
        % Non-dimensional temperature
        Tnd = T(Mask) / Tstar;
        
        % Correlation
        Pnd = 1 - 1.07476 * ( 1 - Tnd.^4.6);
        
        % Allocation and selective assignment
        P = T * 0;
        P(Mask) = Pnd * Pstar;

    else
        P = [];
    end
end