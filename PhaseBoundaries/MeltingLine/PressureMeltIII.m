function P = PressureMeltIII(T)
    
    % Determine the temperatures in the correct range
    Mask = (T > 251.165) & (T < 256.164);
    
    % If the mask is not empty, do the work
    if not(isempty(Mask))
        
        % Nondimensionalizers
        Tstar = 251.165         ; %[K]
        Pstar = 2.085660E+08    ; %[K]
        
        % Non-dimensional temperature
        Tnd = T(Mask) / Tstar;
        
        % Correlation
        Pnd = 1 - 0.299948 * ( 1 - Tnd.^60);
        
        % Allocation and selective assignment
        P = T * 0;
        P(Mask) = Pnd * Pstar;

    else
        P = [];
    end
end