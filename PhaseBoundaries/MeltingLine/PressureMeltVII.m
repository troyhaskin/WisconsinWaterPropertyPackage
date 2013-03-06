function P = PressureMeltVII(T)
    
    % Determine the temperatures in the correct range
    Mask = (T > 355.0) & (T < 715.0);
    
    % If the mask is not empty, do the work
    if not(isempty(Mask))
        
        % Nondimensionalizers
        Tstar = 355.0       ; %[K]
        Pstar = 2.2160E+09  ; %[K]
        
        % Non-dimensional temperature
        Tnd = T(Mask) / Tstar;
        
        % Correlation
        Pnd = exp( 1.73683E+0 * ( 1 - 1./Tnd  )     - ...
                   5.44606E-2 * ( 1 - Tnd.^5  )     + ...
                   8.06106E-8 * ( 1 - Tnd.^22 ))    ;
        
        % Allocation and selective assignment
        P = T * 0;
        P(Mask) = Pnd * Pstar;

    else
        P = [];
    end
end