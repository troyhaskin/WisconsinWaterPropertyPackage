function didTemp = InternalEnergy_Temperature(rho,T,PhaseCheck)
    
    if nargin < 3
        PhaseCheck = true;
    end
    
    Tc     = CriticalTemperature()   ; %[K]
    didtau = InternalEnergy_tau(rho,T,PhaseCheck);
    
    dtaudTemp = -Tc ./ T.^2;
    didTemp   = didtau .* dtaudTemp;
    
end

