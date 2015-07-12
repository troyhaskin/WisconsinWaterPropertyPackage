function delta = DensityRRND(Pnd,tau)

    %   Perform two-phase check
    [Psat,delL,delG] = SaturationStateGivenTauRRND(tau) ;
    isTwoPhi      = abs(Pnd - Psat) < 100*eps()         ;
    isOnePhi      = not(isTwoPhi)                       ;
    aboveCritical = (Psat == 0)                         ;
    delta         = (Psat >  Pnd).*Pnd.*tau        +...
                    (Psat <= Pnd).*delL            ;
    delta         = aboveCritical.*Pnd.*tau             +...
                    (1-aboveCritical).*delta            ;

    %   One-phi solver
    if any(isOnePhi)
        delta(isOnePhi) = NewtonUpdater(...
            @(d,m) updater(d,m,Pnd(isOnePhi),tau(isOnePhi)),delta(isOnePhi),eps(),30)  ;
    end
    
    %   It is two-phase, but quality cannot be determined
    if any(isTwoPhi)
        delta(isTwoPhi) = NaN;
    end
    
    
    function [ddel,rNorm] = updater(delta,mask,Pnd,tau)

        %   Calculate Newton step
        [r,dr] = residual(delta,Pnd(mask),tau(mask))    ;
        ddel  = r./dr                                   ;
        
        %   Guard against negative densities
        isNegative = (delta-ddel)<0;
        while any(isNegative)
            ddel(isNegative) = 0.5*ddel(isNegative);
            isNegative       = (delta-ddel)<0;
        end
        
        
        %   Residual norm
        rNorm = abs(r);
    end
    
    function [r,dr] = residual(delta,Pnd,tau)
        [~,PhiR_d,PhiR_dd] = HelmholtzResidualCombo__d_dd(delta,tau);
        term               = 1 + delta .* PhiR_d                    ;
        r                  = term.*delta - Pnd.*tau                 ;
        dr                 = term + delta.*(PhiR_d + delta.*PhiR_dd);
    end

end





