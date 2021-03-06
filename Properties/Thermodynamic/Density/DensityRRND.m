function delta = DensityRRND(Pnd,tau,x)

    %   Perform two-phase check
    [Psat,delL,delG] = SaturationStateGivenTauRRND(tau) ;
    isTwoPhi      = abs(Pnd - Psat) < 1E-7              ;
    isOnePhi      = not(isTwoPhi)                       ;
    aboveCritical = (Psat == 0)                         ;
    delta         = (Psat >  Pnd).*delG                 +...
                    (Psat <= Pnd).*delL                 ;
    delta         = aboveCritical.*1                    +...
                    (1-aboveCritical).*delta            ;

    %   One-phi solver
    if any(isOnePhi)
        delta(isOnePhi) = NewtonUpdater(...
            @(d,m) updater(d,m,Pnd(isOnePhi),tau(isOnePhi)),delta(isOnePhi),1E-14,30)  ;
    end
    
    %   It is two-phase, but quality cannot be determined
    if any(isTwoPhi)
        
        %   Raw quality
        x               = x(isTwoPhi)                                   ;
        delta(isTwoPhi) = 1./((1-x)./delL(isTwoPhi) + x./delG(isTwoPhi));
        
        %   Borderline liquid saturated
        mask = x < 0;
        if any(mask)
            delta(mask) = delL(mask);
        end
        
        %   Borderline gas saturated
        mask = x > 1;
        if any(mask)
            delta(mask) = delL(mask);
        end
        
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





