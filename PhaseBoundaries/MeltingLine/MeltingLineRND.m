function [Pnd,tau] = MeltingLineRND(n)
    
    if (nargin < 1) || isempty(n)
        n = 100*ones(5,1);
    end
    
    if (numel(n) == 1)
        n = [n;n;n;n;n];
    end
    
    nTotal = cumsum(n);
    
    %   Melting Line transition (double/triple) points
    doubleTaus = CriticalTemperature() ./ [273.16;251.165;256.164;273.31;355.0;715.0];
    
    %   tau vector
    tau = [...
        linspace(doubleTaus(1),doubleTaus(2),n(1))';
        linspace(doubleTaus(2),doubleTaus(3),n(2))';
        linspace(doubleTaus(3),doubleTaus(4),n(3))';
        linspace(doubleTaus(4),doubleTaus(5),n(4))';
        linspace(doubleTaus(5),doubleTaus(6),n(5))'];
    
    pi = tau*0;
    I = 1:n(1);                 pi(I) = PressureMeltIh (doubleTaus(1)./tau(I));
    I = nTotal(1) + (1:n(2));   pi(I) = PressureMeltIII(doubleTaus(2)./tau(I));
    I = nTotal(2) + (1:n(3));   pi(I) = PressureMeltV  (doubleTaus(3)./tau(I));
    I = nTotal(3) + (1:n(4));   pi(I) = PressureMeltVI (doubleTaus(4)./tau(I));
    I = nTotal(4) + (1:n(5));   pi(I) = PressureMeltVII(doubleTaus(5)./tau(I));
    
    Pnd = pi/DimensioningPressure();
    
end

function P = PressureMeltIh(Tnd)
    P = 1 + 1.19539337E6 * ( 1 - Tnd.^3.0000E0) + ...
            8.08183159E4 * ( 1 - Tnd.^2.5750E1) + ...
            3.33826860E3 * ( 1 - Tnd.^1.0375E2) ;
    P = P * 611.657                             ;
end


function P = PressureMeltIII(Tnd)
    P = 1 - 0.299948 * ( 1 - Tnd.^60)   ;
    P = P * 2.085660E+08                ;
end


function P = PressureMeltV(Tnd)
    P = 1 - 1.18721 * ( 1 - Tnd.^8) ;
    P = P * 3.5010E+08              ;
end


function P = PressureMeltVI(Tnd)
    P = 1 - 1.07476 * ( 1 - Tnd.^4.6)   ;
    P = P * 6.3240E+08                  ;
end


function P = PressureMeltVII(Tnd)
        Pnd = exp( 1.73683E+0 * ( 1 - 1./Tnd  ) - ...
                   5.44606E-2 * ( 1 - Tnd.^5  ) + ...
                   8.06106E-8 * ( 1 - Tnd.^22 ));
        P = Pnd * 2.2160E+09                    ;
end
