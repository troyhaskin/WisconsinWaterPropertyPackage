function [Pnd,delL,delG] = SaturationStateGivenTauRRNDGaussSeidel(tau,delL0,delG0,UniqueMask)
    
    [tau,SizeTsat] = Columnify(tau);
    
    % Uniqueness filter; if a UniqueMask is passed; tauSat is assumed to have unique entries.
    if (nargin < 4) || isempty(UniqueMask)
        [tau,~,UniqueMask] = unique(tau);
    end
    
    % Reduced liquid density guess value
    if (nargin < 2) || isempty(delL0)
        delL  = EstimateDelLFromTau(tau);
    else
        delL = delL0;
    end
    
    % Reduced gas density guess value
    if (nargin < 3) || isempty(delG0)
        delG  =  EstimateDelGFromTau(tau);
    else
        delG = delG0;
    end
    
    % Check tau for how close to critical it is
    NearTc     = tau < DoNoIterationTau()   ;
    AboveTc    = tau <   1                  ;   % Check if Tsat is above the critical point
    AboveTt    = tau <= TriplePointTau()    ;
    NotAboveTc = not(AboveTc)               ;
    BelowTt    = not(AboveTt)               ;
    
    
    % Iteration setup
    Calculate = not(NearTc) & NotAboveTc & AboveTt      ; % Mask for temps. not close to the critical point
    
    % Solve the system
    if any(Calculate)
         [delL,delG] = UpdateGaussSeidel(delL,delG,tau);
    end
    
    % These values above for T above the vapor dome ensure that the fluid will
    % be flagged as single phase.
    delL(AboveTc | BelowTt) = 0.0;
    delG(AboveTc | BelowTt) = 0.0;
    
    % If exactly critical, assign exact values
    delL(tau == 1) = 1;
    delG(tau == 1) = 1;
    
    % Saturation pressure
    Pnd             = delG                                              ; % allocate
    Pnd(NotAboveTc) = PressureOneRND(delL(NotAboveTc),tau(NotAboveTc))  ;
    Pnd(tau == 1)   = CriticalPressureND()                              ;
    Pnd(AboveTc | BelowTt) = 0                                          ;
    
    % Expand vectors back to non-unique lengths
    Pnd  = Pnd (UniqueMask) ;
    delL = delL(UniqueMask) ;
    delG = delG(UniqueMask) ;
    
    % Reshape to input shape
    Pnd  = RestoreShape(Pnd,SizeTsat);
    delL = RestoreShape(delL,SizeTsat);
    delG = RestoreShape(delG,SizeTsat);
    
end




function [delL,delG] = UpdateGaussSeidel(delL,delG,tau)
    
    Tolerance = 1E-6;
    
    PhiR    = @(delta,M) HelmholtzResidual   (delta,tau(M));
    PhiR_d  = @(delta,M) HelmholtzResidual_d (delta,tau(M));
    
    NotConverged = true();
    
    Psig = @(delL,delG,M) (PhiR(delL,M) - PhiR(delG,M) + log(delL./delG)).*delL.*delG ;
    RL   = @(delL,delG,M) delL .* (delL-delG) .* (1 + delL .* PhiR_d(delL,M)) - Psig(delL,delG,M);
    RG   = @(delL,delG,M) delG .* (delL-delG) .* (1 + delG .* PhiR_d(delG,M)) - Psig(delL,delG,M);
    dRL  = @(delL,delG,M) (RL(delL + 1E-8 , delG,M) - RL(delL - 1E-8 , delG,M))./(2E-8);
    dRG  = @(delL,delG,M) (RG(delL , delG + 1E-8,M) - RG(delL , delG - 1E-8,M))./(2E-8);
    
    All  = true(length(tau),1);
    Rk   = RL(delL,delG,All);
    Mask = abs(Rk) > Tolerance;
    while NotConverged
        
        while any(Mask)
            delL(Mask) = delL(Mask) - Rk(Mask)./dRL(delL(Mask),delG(Mask),Mask);
            Rk(Mask)   = RL(delL(Mask),delG(Mask),Mask);
            Mask = abs(Rk) > Tolerance;
        end
        Rk = RG(delL,delG,All);
        Mask = abs(Rk) > Tolerance;
        while any(Mask)
            delG(Mask) = delG(Mask) - Rk(Mask)./dRG(delL(Mask),delG(Mask),Mask);
            Rk(Mask) = RG(delL(Mask),delG(Mask),Mask);
            Mask = abs(Rk) > Tolerance;
        end

        Rk = RL(delL,delG,All);
        Mask = abs(Rk) > Tolerance;
        NotConverged = any(abs(Rk) > Tolerance);

        
        Show(Rk,'%+10.6E');
        
    end
    
end