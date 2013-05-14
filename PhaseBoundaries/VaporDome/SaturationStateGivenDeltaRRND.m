function [Pnd,tau,delL,delG] = SaturationStateGivenDeltaRRND(delGiven,UniqueMask)
    
    % ========================================================= %
    %                        Set-Up                             %
    % ========================================================= %
    
    % Filter out non-unique entries
    if (nargin < 2) || isempty(UniqueMask)
        [delGiven,~,UniqueMask] = unique(delGiven,'stable');
    end

    % Given Masks
    GivenL = delGiven >= 1;
    GivenG = delGiven <= 1;

    % Ordering Masks: since the solution algorithm relies on the ordering of dels to be
    % [delL;delG], these index arrays allow insertion of the solved properties into the 
    % original ordering.
    Ioriginal = 1 : length(delGiven);
    Iordered  = [Ioriginal(GivenL),Ioriginal(GivenG)];

    % Given densities
    delLgiven = delGiven(GivenL)        ;
    delGgiven = delGiven(GivenG)        ;
    delGiven  = [delLgiven;delGgiven]   ; % destroys original ordering; use index arrays above

    % Number of gas densities given
    NliquidGiven = length(delLgiven);



    % ================================================================================== %
    %                                Computability Filters                               %
    % ================================================================================== %
    %   There are special considerations for solving the saturation system with del as the
    %   independent parameter:
    %       o   There is a limit on the minimum value of liquid density and maximum value 
    %           of gas density accepted for iteration.  These limits occur near the critical
    %           point where the saturation system becomes increasingly ill-conditioned and,
    %           therefore, no iteration will currently be attemped on values falling into 
    %           the forbidden region. This is also done in the tau saturation function.
    %           The values of these limits are retrieved from NearCriticalStabilityDensityLimitsR().
    %
    %       o   There is a maximum liquid density above which and a minimum gas density below
    %           which saturation cannot occur.  The maximum liquid density occurs around
    %           277.15 K, and the minimum gas density occurs on the triple line.
    %           The values of these limits are retrieved from SaturationLineDensityBoundsR().
    %
    %       o   Near the triple point, the saturation line becomes a multi-valued function
    %           of density; this causes numerous calculation issues, both in iteration and
    %           guess value functions.  Current implementation to deal with this is:
    %               1.  Use curve fits to find the average tau at which the del occurs 
    %               2.  Return the saturation values for that tau value
    %               3.  Throw a warning that the function is multi-valued, an average value
    %                   has been returned, and using a function designed to use another 
    %                   property to uniquely specify the saturation state is required.
    %
    
    [delLmax   , delGmin]   = SaturationLineDensityBoundsR()        ;
    [delLnearC , delGnearC] = NearCriticalStabilityDensityLimitsR() ;
    [delLnearT , delGnearT] = NearTripleStabilityDensityLimitsR()   ;

    CanSaturate    = [delLgiven <= delLmax     ;...   % Maximum saturation density bound
                      delGgiven >= delGmin    ];      % Minimum saturation density bound
    DoNotIterate   = [delLgiven <= delLnearC   ;...   % Near critical point stability bound
                      delGgiven >= delGnearC  ];      % Near critical point stability bound
    NotNearTriple  = [delLgiven <= delLnearT   ;...   % Near triple   point stability bound
                      delGgiven >= delGnearT  ];      % Near triple   line  stability bound
    NearTriple     = not(NotNearTriple)        ;
    CannotSaturate = not(CanSaturate)          ;

    % Logical index for values that will be iterated upon
    WillGuess    = NotNearTriple & CanSaturate          ;
    WillIterate  = WillGuess     & not(DoNotIterate)    ;



    % ================================================================================== %
    %                                    Guess Values                                    %
    % ================================================================================== %

    % Get guess values for all given density values not near the triple line
    if any(WillGuess)
        [delLguess,delGguess,tauGuess] = EstimateDelLDelGTauFromDel(delGiven(WillGuess));
    end

    % Get average tau value for near triple point liquid densities
    if any(NearTriple)
        [tauHigh,tauLow]   = GetTauLimitsNearTriplePoint(delGiven(NearTriple))      ;
        tauAverage         = (tauHigh + tauLow)/2                                   ;
        [~,delLnearT,delGnearT] = SaturationStateGivenTauRRND(tauAverage)           ;
        
        warning('WWPP:SaturationStateGivenDeltaRRND:MultiValuedRegime'                  ,...
                ['Some of the supplied densities supplied fall into a regime where '    ,...
                 'the saturation state is not uniquely defined.  Returned values are '  ,...
                 'come the average of the two temperatures at which the density '       ,...
                 'occurs.  If this is not acceptable, please use a function that '      ,...
                 'uses another indepedent property to fully define the saturation '     ,...
                 'state.']);
    else
        delLnearT  = [];
        delGnearT  = [];
        tauAverage = [];
    end



    % ================================================================================== %
    %                                 Newton Solution                                    %
    % ================================================================================== %
    
    if any(WillIterate)

        % Iteration parameters
        Tolerance   = 1E-10;
        IterMax     = 500;

        % Full guess matrix
        Guess = [[delGguess(1:NliquidGiven);delLguess(NliquidGiven+1:end)],...
                   tauGuess(WillIterate)];

        % Update handle
        Updater  = @(Unknowns,Mask) UpdateSystem(Unknowns,Mask,delGiven,NliquidGiven);

        % Newton solution
        Solution = NewtonUpdater(Updater,Guess,Tolerance,IterMax,true);

        % Pull solved dels and taus
        delGsol = Solution(WillIterate & GivenL,1);
        delLsol = Solution(WillIterate & GivenG,1);
        tauSol  = Solution(:,2);

    else
        % Initialize as empty
        delGsol = [];
        delLsol = [];
        tauSol  = [];
    end
    
    % Allocations
    delG = delGiven;
    delL = delGiven;
    tau  = delGiven;
    
    % Insert guess values
    delG(WillGuess & GivenL) = delGguess  ;
    delL(WillGuess & GivenG) = delLguess  ;
    tau (WillGuess)          = tauGuess   ;
    
    % Insert Newton solutions
    delG(WillIterate & GivenL) = delGsol    ;
    delL(WillIterate & GivenG) = delLsol    ;
    tau (WillIterate)          = tauSol     ;

    % Insert Near Triple results
    delG(NearTriple) = delGnearT    ;
    delL(NearTriple) = delLnearT    ;
    tau (NearTriple) = tauAverage   ;

    % Undefinable quantities
    delG(CannotSaturate) = NaN  ;
    delL(CannotSaturate) = NaN  ;
    tau (CannotSaturate) = NaN  ;

    % Put back into passed-in ordering
    delG(Iordered) = delG;
    delL(Iordered) = delL;
    tau (Iordered) = tau;

    % Get Pressures
    Pnd = PressureOneRND(delL,tau);

    % Splay out uniques
    delG = delG(UniqueMask);
    delL = delL(UniqueMask);
    tau  = tau (UniqueMask);
    Pnd  = Pnd (UniqueMask);

end





function [dUnknowns,RNorm] = UpdateSystem(Unknowns,Mask,delGiven,NliquidGiven)
    
    % Pull unknowns
    dels = Unknowns(:,1);
    taus = Unknowns(:,2);
    
    % Logical masks for delGiven
    lGivenL = (Mask <= NliquidGiven)  ; % liquid density given
    lGivenG = not(lGivenL)            ; % gas    density given
    
    % Integer masks for the Unknowns:
    %       These are different from the delGiven masks since the Unknowns are returned 
    %       to this function already contracted (i.e., masked).
    %
    Nmask   = length(Mask)          ; % Number of unconverged values
    NgivenL = sum(lGivenL)          ; % number of liquid knowns remaining
    iGivenL =       1    : NgivenL  ; % integer mask for unknown delGs
    iGivenG = (NgivenL+1): Nmask    ; % integer mask for unknown delLs

    % Assign the Unknowns to descriptive variables
    delGgivenL = dels(iGivenL);
    delLgivenG = dels(iGivenG);
    tauLgivenL = taus(iGivenL);
    tauGgivenG = taus(iGivenG);

    % Pull delGiven for unconverged Unknowns
    delLgivenL = delGiven(lGivenL);
    delGgivenG = delGiven(lGivenG);

    % Combined vectors for liquid/gas-given agnostic (LGA) evaluations (aides in vectorization)
    delL = [delLgivenL ; delLgivenG ];
    delG = [delGgivenL ; delGgivenG ];
    tau  = [tauLgivenL ; tauGgivenG ];

    % Integer masks for vectorized Helmholtz functions:
    %       The most expensive evaluation in almost all of these thermodynamic calls is 
    %       the Helmholtz free energy functions.  As such, one the biggest optimizations
    %       lies in reducing calls to those functions by vectorizing as much as possible.
    %       


    % Helmholtz energy function handles (d = del, t = tau; shortened to not 
    % conflict with above variables).
    PhiR    = @(d,t) HelmholtzResidual   (d,t);
    PhiR_d  = @(d,t) HelmholtzResidual_d (d,t);
    PhiR_t  = @(d,t) HelmholtzResidual_t (d,t);
    PhiR_dd = @(d,t) HelmholtzResidual_dd(d,t);
    PhiR_dt = @(d,t) HelmholtzResidual_dt(d,t);


    % ============================================================================ %
    %                      Calculate Dimensionless Pressure                        %
    % ============================================================================ %
    % Pnd is the dimensionless pressure eliminated from
    % the system to shrink the solution space.  Here are the functions and
    % derivatives that are needed in the residuals below.

    % Pnd = Pnd1 * Pnd2 to more easily maintain/check the handles.
    Pnd1 = PhiR(delL,tau) - PhiR(delG,tau) + log(delL./delG); % LGA
    Pnd2 = (delL.*delG)./(delL - delG)                      ; % LGA

    % tau derivatives  (LGA)
    Pnd1_t = PhiR_t(delL,tau) - PhiR_t(delG,tau);
    Pnd2_t = 0                                  ;

    % del derivatives while solving for delG, tauG with delL given
    if any(GivenL)
        Pnd1_dL = 1./delLL + PhiR_d(delLL,tauL) ;
        f       = delGL./(delLL - delGL)        ;
        Pnd2_dL = (1 - delLL./(delLL-delGL)).*f ;
    else
        Pnd1_dL = [];
        Pnd2_dL = [];
    end
    
    % del derivatives while solving for delL, tauL with delG given
    if any(GivenG)
        Pnd1_dG = -1./delGG - PhiR_d(delGG,tauG);
        f        = delLG./(delLG - delGG)       ;
        Pnd2_dG = (1 + delLG./(delLG-delGG)).*f ;
    else
        Pnd1_dG = [];
        Pnd2_dG = [];
    end
    
    % Create combined Pnd*_d vectors
    Pnd1_d = [Pnd1_dL ; Pnd1_dG];
    Pnd2_d = [Pnd2_dL ; Pnd2_dG];
    
    % Pnd calculated
    Pnd    = Pnd1   .* Pnd2                    ;
    Pnd_d  = Pnd1_d .* Pnd2 + Pnd2_d .* Pnd1   ;
    Pnd_t  = Pnd1_t .* Pnd2 + Pnd2_t .* Pnd1   ;

    % Residuals
    R1   = Pnd  - delL .* (1 + delL.*PhiR_d(delL,tau));
    R2   = Pnd  - delG .* (1 + delG.*PhiR_d(delG,tau));
    
    % Jacobian values
    R1_d = Pnd_d - (1 + 2*delL.*PhiR_d(delL,tau) + delL.^2.*PhiR_dd(delL,tau))  ;
    R1_t = Pnd_t - delL.^2 .* PhiR_dt(delL,tau)                                 ;
    R2_d = Pnd_d                                                                ;
    R2_t = Pnd_t - delG.^2 .* PhiR_dt(delG,tau)                                 ;
    DetR = R1_d .* R2_t - R1_t.*R2_d                                            ;

    % Final Newton Directions
    ddel  =  (R1 .* R2_t - R2 .* R1_t) ./ DetR ;
    dtau  =  (R2 .* R1_d - R1 .* R2_d) ./ DetR ;

    % Final outputs
    dUnknowns = [ddel,dtau]         ;
    RNorm     = abs(R1) + abs(R2)   ; % L_1 norm

end

function [PhiRgivenL,PhiRgivenG] = HelmholtzHarvest(delL,delG,tauL,tauG,GivenL,GivenG)
    
    Nliquid = length(delL);
    dels    = [delL;delG];
    taus    = [tauL;tauG];
    PhiR    = HelmholtzResidual(dels,taus);
    
    
    PhiRgivenL.L
    PhiRgivenL.L_d
    PhiRgivenL.L_t
    PhiRgivenL.L_dt
    PhiRgivenL.L_dd
    PhiRgivenL.G
    PhiRgivenL.G_d
    PhiRgivenL.G_t
    PhiRgivenL.G_dt
    PhiRgivenL.G_dd

    PhiRgivenG.G
    PhiRgivenG.G_d
    PhiRgivenG.G_t
    PhiRgivenG.G_dt
    PhiRgivenG.G_dd
    PhiRgivenG.L
    PhiRgivenG.L_d
    PhiRgivenG.L_t
    PhiRgivenG.L_dt
    PhiRgivenG.L_dd
end



