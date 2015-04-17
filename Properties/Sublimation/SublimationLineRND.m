function [Pnd,tau] = SublimationLineRND(n)
    
    if (nargin < 1) || isempty(n)
        n = 100;
    end
    
    %   Temperature vectors
    Tnd = linspace(50,273.16,n)'    ;
    tau = CriticalTemperature()./Tnd;
    Tnd = Tnd/273.16                ;
    
    % Correlation for dimensionless pressure
    pi = exp((...
        0.273203819E+2 * Tnd.^0.120666667E+1  -  ...
        0.212144006E+2 * Tnd.^0.333333333E-2  -  ...
        0.610598130E+1 * Tnd.^0.170333333E+1) ./ Tnd);
    
    %   Consistently non-dimensionalize
    Pnd = 611.657 * pi / DimensioningPressure();

end