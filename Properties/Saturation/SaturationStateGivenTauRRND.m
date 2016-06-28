function [Pnd,delL,delG] = SaturationStateGivenTauRRND(tau,delL0,delG0,UniqueMask)
    
    % Uniqueness filter; if a UniqueMask is passed; tauSat is assumed to have unique entries.
    if (nargin < 4) || isempty(UniqueMask)
        [tau,~,UniqueMask] = unique(tau);
    end
    
    % Reduced liquid density guess value
    if (nargin < 2) || isempty(delL0) || any(delL0 <= 0)
        delL  = EstimateDelLFromTau(tau);
    else
        delL = delL0;
    end
    
    % Reduced gas density guess value
    if (nargin < 3) || isempty(delG0) || any(delG0 <= 0)
        delG  =  EstimateDelGFromTau(tau);
    else
        delG = delG0;
    end
    
    % Check tau for how close to critical it is
    nearTc     = tau < DoNoIterationTau()   ;
    aboveTc    = tau <   1                  ;   % Check if Tsat is above the critical point
    aboveTt    = tau <= TriplePointTau()    ;
    notAboveTc = not(aboveTc)               ;
    belowTt    = not(aboveTt)               ;
    
    
    % Iteration setup
    Tolerance = DefaultAbsoluteIterationTolerance()     ; % Abolsute iteration tolerance
    IterMax   = DefaultMaximumIterationCount()          ; % Maximum iteration count
    Calculate = not(nearTc) & notAboveTc & aboveTt      ; % Mask for temps. not close to the critical point
    
    % Solve the system
    if any(Calculate)
        
        %   Perform one Newton calculation to initialize Broyden Jacobian
        Guess     = [delL(Calculate),delG(Calculate)]       ; % Starting values for the iteration
        [dx,~,S.r1,S.r2,S.iJ11,S.iJ12,S.iJ21,S.iJ22] = newton(Guess,1:nnz(Calculate),tau(Calculate));
        S.delL    = Guess(:,1);
        S.delG    = Guess(:,2);
        
        %   Solve
%         xSol = NewtonUpdater(@(x,mask) broydenClose(x,mask,tau(Calculate)),Guess-dx,1E-13,IterMax);
        xSol = NewtonUpdater(@(x,mask)       newton(x,mask,tau(Calculate)),Guess,Tolerance,IterMax);
        
        % Update the iterated values
        delL(Calculate) = xSol(:,1);
        delG(Calculate) = xSol(:,2);
        
    end
    
    % These values above for T above the vapor dome ensure that the fluid will
    % be flagged as single phase.
    delL(aboveTc | belowTt) = 0.0;
    delG(aboveTc | belowTt) = 0.0;
    
    % If exactly critical, assign exact values
    delL(tau == 1) = 1;
    delG(tau == 1) = 1;
    
    % Saturation pressure
    Pnd                    = delG                                               ;
    Pnd(notAboveTc)        = PressureOneRND(delG(notAboveTc),tau(notAboveTc))   ;
    Pnd(tau == 1)          = CriticalPressureND()                               ;
    Pnd(aboveTc | belowTt) = 0                                                  ;
    
    % Expand vectors back to non-unique lengths
    Pnd  = Pnd (UniqueMask) ;
    delL = delL(UniqueMask) ;
    delG = delG(UniqueMask) ;

    
    
    
    function [dx,RNorm] = broydenClose(x,mask,tau0)
        [dx,RNorm,S] = broyden(x,mask,tau0,S);
    end
    
    
end

function [dx,RNorm,r1,r2,iJ11,iJ12,iJ21,iJ22] = newton(x,Mask,tau0)
    
    N    = length(Mask);
    tau  = [tau0(Mask);tau0(Mask)];
    delL = x(:,1);
    delG = x(:,2);
    
    
    %   Helmholtz Free Energy values
    [PhiR,PhiR_d,PhiR_dd] = HelmholtzResidualCombo__d_dd([delL;delG],tau);
    [PhiRL,PhiRG,PhiR_dL,PhiR_dG,PhiR_ddL,PhiR_ddG] = VectorChunk([PhiR;PhiR_d;PhiR_dd],N);
    
    %   Terms used to form the residuals
    Psig1    =  PhiRL - PhiRG + log(delL./delG) ;
    Psig1_dL =   PhiR_dL + 1./delL              ;
    Psig1_dG = -(PhiR_dG + 1./delG)             ;
    
    PdelL = 1 + delL.*PhiR_dL;
    PdelG = 1 + delG.*PhiR_dG;
    
    PdelL_d = PhiR_dL + delL.*PhiR_ddL;
    PdelG_d = PhiR_dG + delG.*PhiR_ddG;
    
    delLmG = (delL-delG);
    
    
    %   Residual and partial derivative values
    r1 = delG.*Psig1 - PdelL.*delLmG;
    r2 = delL.*Psig1 - PdelG.*delLmG;
    
    r1_dL = delG.*Psig1_dL - (delLmG.*PdelL_d + PdelL)  ;
    r1_dG = Psig1 + delG.*Psig1_dG + PdelL              ;
    
    r2_dL = Psig1 - PdelG  + delL.*Psig1_dL;
    r2_dG = delL.*Psig1_dG - (delLmG.*PdelG_d - PdelG);

    
    % Determinant
    DetJ = r1_dL .* r2_dG - r1_dG .*r2_dL;
    
    % Inverse Jacobian elements
    iJ11 =  r2_dG ./ DetJ;
    iJ12 = -r1_dG ./ DetJ;
    iJ21 = -r2_dL ./ DetJ;
    iJ22 =  r1_dL ./ DetJ;


    % Newton updates
    ddelL = iJ11 .* r1 + iJ12 .* r2;
    ddelG = iJ21 .* r1 + iJ22 .* r2;
    
    % Pack updates and calculate norm for the Newton updater
    dx    = [ddelL,ddelG]       ;
    RNorm = abs(r1) + abs(r2)   ;

end


function [dx,RNorm,S] = broyden(xk,mask,tau0,S)
    
    N    = length(mask);
    tau  = [tau0(mask);tau0(mask)];
    delL = xk(:,1);
    delG = xk(:,2);
    
    
    %   Helmholtz Free Energy values
    [PhiR,PhiR_d]                 = HelmholtzResidualCombo__d([delL;delG],tau)  ;
    [PhiRL,PhiRG,PhiR_dL,PhiR_dG] = VectorChunk([PhiR;PhiR_d],N)                ;
    
    %   Terms used to form the residuals
    Psig1  = PhiRL - PhiRG + log(delL./delG);
    PdelL  = 1 + delL.*PhiR_dL              ;
    PdelG  = 1 + delG.*PhiR_dG              ;

    %   Residual and partial derivative values
    r1 = delG.*Psig1 - PdelL.*(delL-delG);
    r2 = delL.*Psig1 - PdelG.*(delL-delG);
    
    % Form deltas
    dx1  = delL - S.delL(mask)  ;
    dx2  = delG - S.delG(mask)  ;
    dr1  = r1   - S.r1(mask)    ;
    dr2  = r2   - S.r2(mask)    ;
    drN  = dr1.^2 + dr2.^2      ;
    
    % Update Inverse Jacobian
    term         = (dx1 - S.iJ11(mask).*dr1 - S.iJ12(mask).*dr2)./drN;
    S.iJ11(mask) = S.iJ11(mask) + dr1.*term;
    S.iJ12(mask) = S.iJ12(mask) + dr2.*term;
    term         = (dx2 - S.iJ21(mask).*dr1 - S.iJ22(mask).*dr2)./drN;
    S.iJ21(mask) = S.iJ21(mask) + dr1.*term;
    S.iJ22(mask) = S.iJ22(mask) + dr2.*term;


    % Newton updates
    ddelL = S.iJ11(mask) .* r1 + S.iJ12(mask) .* r2;
    ddelG = S.iJ21(mask) .* r1 + S.iJ22(mask) .* r2;
    
    % Pack updates and calculate norm for the Newton updater
    dx    = [ddelL,ddelG]       ;
    RNorm = abs(r1) + abs(r2)   ;
    
    
    %   Update struct elements
    S.delL(mask) = delL     ;
    S.delG(mask) = delG     ;
    S.r1(mask)   = r1       ;
    S.r2(mask)   = r2       ;

end

