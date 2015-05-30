function [Deltabi_d,Deltabi_t,Deltabi_dt,Deltabi_dd] = ...
        GetDeltabiCombo_d_t_dt_dd(delta,deltaMod,Delta,Theta,A,B,a,b,betaInv)
    
    %   Deltas and Parts
    Part        = 4 * A .* Theta .* betaInv .* (0.5*betaInv-1);
    Part        = Part .* deltaMod.^(0.5*betaInv-2);
    Part        = Part + 4 * B .* a .* (a-1) .* deltaMod.^(a-2);
    Part        = Part + 2 * A.^2 .* betaInv.^2 .* deltaMod.^(betaInv-2);
    Part        = deltaMod .* Part;
    Delta_d     = GetDelta_d (delta,deltaMod,Theta,A,B,a,betaInv);
    Delta_dd    = Part + Delta_d./(delta-1);
    Partt       = A * betaInv.*Delta.*(delta-1).* deltaMod.^(0.5*betaInv-1);
    
    
    %   Deltabis
    Deltabi_d   = b    * Delta.^(b-1).*Delta_d;
    Deltabi_t   = -2*b * Theta .* Delta.^(b-1);
    Deltabi_dt  = -2*b*Delta.^(b-2).* (Partt + (b-1)*Theta.*Delta_d);
    Deltabi_dd  = b    * (Delta.^(b-1).*Delta_dd + (b-1)*Delta.^(b-2).*Delta_d.^2) ;
end
