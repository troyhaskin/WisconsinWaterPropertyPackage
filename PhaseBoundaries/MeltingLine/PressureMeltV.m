function P = PressureMeltV(T)
    
    % Determine the temperatures in the correct range
    Mask = (T > 256.164) & (T < 273.31);
    
    % If the mask is not empty, do the work
    if not(isempty(Mask))
        
        % Nondimensionalizers
        Tstar = 256.164     ; %[K]
        Pstar = 3.5010E+08  ; %[K]
        
        % Non-dimensional temperature
        Tnd = T(Mask) / Tstar;
        
        % Correlation
        Pnd = 1 - 1.18721 * ( 1 - Tnd.^8);
        
        % Allocation and selective assignment
        P = T * 0;
        P(Mask) = Pnd * Pstar;

    else
        P = [];
    end
end