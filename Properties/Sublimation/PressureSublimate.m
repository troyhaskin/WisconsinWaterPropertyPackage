function P = PressureSublimate(T)
    
    % Determine the temperatures in the correct range
    Mask = (T >= 50.0) & (T <= 273.16);
    
    % If the mask is not empty, do the work
    if not(isempty(Mask))
        
        % Nondimensionalizers
        Tstar = 273.16      ; %[K]
        Pstar = 611.657     ; %[Pa]
        
        % Non-dimensional temperature
        Tnd = T(Mask) / Tstar;
        
        % Correlation
        Pnd = exp( (2.73203819E+1 * Tnd.^1.20666667E+0  -  ...
                    2.12144006E+1 * Tnd.^3.33333333E-3  -  ...
                    6.10598130E+0 * Tnd.^1.70333333E+0)  ./ ...
                    Tnd)    ;
        
        % Allocation and selective assignment
        P = T * 0;
        P(Mask) = Pnd * Pstar;

    else
        P = [];
    end
end