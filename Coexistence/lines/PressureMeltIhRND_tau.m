function Pnd_tau = PressureMeltIhRND_tau(tau)
    
    %   Convert to fit's non-dimensional temperature
    Tnd = TriplePointTau()./tau;
    
    %   Derivative of correlation for dimensionless pressure w.r.t. Tnd
    pi_Tnd =    1.19539337E6 * -3.0000E+0 * Tnd.^2.0000E0   + ...
                8.08183159E4 * -2.5750E+1 * Tnd.^2.4750E1   + ...
                3.33826860E3 * -1.0375E+2 * Tnd.^1.0275E2   ;
    
    %   Convert to WWPP non-dimensional pressure
    Tnd_tau = -Tnd./tau;
    Pnd_tau = 611.657/DimensioningPressure() * pi_Tnd .* Tnd_tau;
    
end