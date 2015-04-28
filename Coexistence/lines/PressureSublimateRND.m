function [Pnd,Pnd_tau] = PressureSublimateRND(tau)
    
    %   Convert to fit's non-dimensional temperature
    Tnd = TriplePointTau()./tau;
    
    %   Correlation for dimensionless pressure
    pi = exp((...
        0.273203819E+2 * Tnd.^0.120666667E+1  -  ...
        0.212144006E+2 * Tnd.^0.333333333E-2  -  ...
        0.610598130E+1 * Tnd.^0.170333333E+1) ./ Tnd);
    
    %   Convert to WWPP non-dimensional pressure
    Pnd = 611.657/DimensioningPressure() * pi;


    % Return derivative if requested
    if (nargout > 1)
        pi_Tnd = pi .* (...
            0.273203819E+2 * (0.120666667E+1 - 1) * Tnd.^0.120666667E+1  -  ...
            0.212144006E+2 * (0.333333333E-2 - 1) * Tnd.^0.333333333E-2  -  ...
            0.610598130E+1 * (0.170333333E+1 - 1) * Tnd.^0.170333333E+1) ./ Tnd^2;
    
        %   Convert to WWPP non-dimensional pressure
        Tnd_tau = -Tnd./tau;
        Pnd_tau = 611.657/DimensioningPressure() * pi_Tnd * Tnd_tau;
    end

end