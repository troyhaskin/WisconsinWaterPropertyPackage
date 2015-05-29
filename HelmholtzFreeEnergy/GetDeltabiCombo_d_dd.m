function [Deltabi_d,Deltabi_dd] = GetDeltabiCombo_d_dd(delta,deltaMod,Delta,Theta,A,B,a,b,betaInv)
    
    %   Deltas
    Part        = 4 * A .* Theta .* betaInv .* (0.5*betaInv-1);
    Part        = Part .* deltaMod.^(0.5*betaInv-2);
    Part        = Part + 4 * B .* a .* (a-1) .* deltaMod.^(a-2);
    Part        = Part + 2 * A.^2 .* betaInv.^2 .* deltaMod.^(betaInv-2);
    Part        = deltaMod .* Part;
    Delta_d     = GetDelta_d (delta,deltaMod,Theta,A,B,a,betaInv);
    Delta_dd    = Part + Delta_d./(delta-1);
    
    %   Deltabis
    Deltabi_d   = b * Delta.^(b-1).*Delta_d;
    Deltabi_dd  = b * (Delta.^(b-1).*Delta_dd + (b-1)*Delta.^(b-2).*Delta_d.^2) ;
end
