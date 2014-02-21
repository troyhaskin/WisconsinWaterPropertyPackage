function [Pnd,delL,delG] = SaturationStateGivenTauRRND(tau,delL0,delG0,UniqueMask)
    
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
    Tolerance = 1E-12                                   ; % Abolsute iteration tolerance
    IterMax   = 1E3                                     ; % Maximum iteration count
    Calculate = not(NearTc) & NotAboveTc & AboveTt      ; % Mask for temps. not close to the critical point
    Guess     = [delL(Calculate),delG(Calculate)]       ; % Starting values for the iteration
    UpdateFun = @(x,Mask) Updater(x,Mask,tau(Calculate)); % Function used for updating the solution
    
    % Solve the system
    if any(Calculate)
        xSol = NewtonUpdater(UpdateFun,Guess,Tolerance,IterMax);

        % Update the iterated values
        delL(Calculate) = xSol(:,1);
        delG(Calculate) = xSol(:,2);

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


function [dx,RNorm] = Updater(x,Mask,tau0)
    
    tau  = tau0(Mask);
    delL = x(:,1);
    delG = x(:,2);

    PhiR    = @(delta) HelmholtzResidual   (delta,tau);
    PhiR_d  = @(delta) HelmholtzResidual_d (delta,tau);
    PhiR_dd = @(delta) HelmholtzResidual_dd(delta,tau);
    
    Psig1    = PhiR(delL) - PhiR(delG) + log(delL./delG)    ;
    Psig1_dl =  PhiR_d(delL) + 1./delL                      ;
    Psig1_dg =-(PhiR_d(delG) + 1./delG)                     ;
    
    PdelL = 1 + delL.*PhiR_d(delL);
    PdelG = 1 + delG.*PhiR_d(delG);
    
    PdelL_d = PhiR_d(delL) + delL.*PhiR_dd(delL);
    PdelG_d = PhiR_d(delG) + delG.*PhiR_dd(delG);
    
    R1 = delG.*Psig1 - PdelL.*(delL-delG);
    R2 = delL.*Psig1 - PdelG.*(delL-delG);
    
    R1_dl = delG.*Psig1_dl + (delG-delL).*PdelL_d - PdelL;
    R1_dg = Psig1 + PdelL + delG.*Psig1_dg;
    
    R2_dl = Psig1 - PdelG + delL.*Psig1_dl;
    R2_dg = delL.*Psig1_dg + (delG-delL).*PdelG_d + PdelG;
    
    DetJ = R1_dl .* R2_dg - R1_dg .*R2_dl;
    
    iJ11 =  R2_dg ./ DetJ;
    iJ12 = -R1_dg ./ DetJ;
    iJ21 = -R2_dl ./ DetJ;
    iJ22 =  R1_dl ./ DetJ;
    
    ddelL = iJ11 .* R1 + iJ12 .* R2;
    ddelG = iJ21 .* R1 + iJ22 .* R2;
    
    dx    = [ddelL,ddelG];
    RNorm = abs(R1) + abs(R2);
    
end
