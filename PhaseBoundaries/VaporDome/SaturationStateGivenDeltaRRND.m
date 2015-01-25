function [Pnd,tau,delL,delG] = SaturationStateGivenDeltaRRND(delGiven,UniqueMask)
    
    % ========================================================= %
    %                        Set-Up                             %
    % ========================================================= %
    
    % Filter out non-unique entries
    if (nargin < 2) || isempty(UniqueMask)
        [delGiven,~,UniqueMask] = UniqueEnough(delGiven,min(eps(delGiven)));
    end

    % Given Masks
    GivenL = delGiven >= 1;
    GivenG = delGiven <= 1;

    % Ordering Masks: since the solution algorithm relies on the ordering of dels to be
    % [delL;delG], these index arrays allow insertion of the solved properties into the 
    % original ordering.
    Ngiven    = length(delGiven);
    Ioriginal = 1 : Ngiven;
    Iordered  = [Ioriginal(GivenL),Ioriginal(GivenG)];

    % Given densities
    delLgiven = delGiven(GivenL)        ;
    delGgiven = delGiven(GivenG)        ;
    delGiven  = [delLgiven;delGgiven]   ; % destroys original ordering; use index arrays above

    % Number of gas densities given
    NgivenL = length(delLgiven) ;
    NgivenG = Ngiven - NgivenL  ;
    lGivenL = [ true(NgivenL,1);false(NgivenG,1)];
    lGivenG = [false(NgivenL,1); true(NgivenG,1)];


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
    
    [delLmax   , delGmin  ] = SaturationLineDensityBoundsR()        ;
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
                ['Some of the supplied densities fall into a regime where '    ,...
                 'the saturation state is not uniquely defined.  Returned values are '  ,...
                 'the average of the two temperatures at which the density '            ,...
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
        Tolerance   = DefaultAbsoluteIterationTolerance();
        IterMax     = DefaultMaximumIterationCount();

        % Fill guess matrix
        Guess = [[delGguess(1:NgivenL);delLguess(NgivenL+1:end)],...
                   tauGuess(WillIterate),delGiven];

        % Update handle
        Updater  = @(Unknowns,Mask) UpdateSystem(Unknowns,Mask,NgivenL);

        % Newton solution
        Solution = NewtonUpdater(Updater,Guess,Tolerance,IterMax,true);

    end
    
    % Allocations
    delG = delGiven;
    delL = delGiven;
    tau  = delGiven;
    
    % Insert non-iterated guess values
    delG(not(WillIterate) & lGivenL) = delGguess(not(WillIterate))   ;
    delL(not(WillIterate) & lGivenG) = delLguess(not(WillIterate))   ;
    tau (not(WillIterate))           = tauGuess (not(WillIterate))   ;
    
    % Insert Newton solutions
    delL(WillIterate & lGivenG) = Solution(lGivenG,1)   ;
    delG(WillIterate & lGivenL) = Solution(lGivenL,1)   ;
    tau (WillIterate)           = Solution(   :   ,2)   ;

    % Insert Near Triple results
    delG(NearTriple) = delGnearT    ;
    delL(NearTriple) = delLnearT    ;
    tau (NearTriple) = tauAverage   ;

    % Undefinable quantities
    delG(CannotSaturate) = 0  ;
    delL(CannotSaturate) = 0  ;
    tau (CannotSaturate) = 0  ;

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





function [dUnknowns,RNorm] = UpdateSystem(Unknowns,Mask,NliquidGiven)
%
%   This function provides the Newton updates for an unknown reduced density
%   delUnknown and associated saturation temperature tauSat given the opposing reduced 
%   density delKnown.  The residual to be minimized is
%       --                  --   --                                                  --
%       | RL(delL,delG,tauS) |   | delL [f delG - (delL - delG) (1 + delL PhiR_delL)] |
%       |                    | = |                                                    |
%       | RG(delL,delG,tauS) |   | delG [f delL - (delL - delG) (1 + delG PhiR_delG)] |
%       --                  --   --                                                  --
%   where
%
%       f         = PhiR(delL,tauS) - PhiR(delG,tauS) + Log(delL/delG)
%       PhiR_del* = PhiR_del(del*,tauS)
%
%
    
    % Pull unknowns
    dels     = Unknowns(:,1)    ;
    taus     = Unknowns(:,2)    ;
    delGiven = Unknowns(:,3)    ;
    
    Nmask   = length(Mask)              ;
    NgivenL = nnz(Mask <= NliquidGiven) ; % number of liquid knowns remaining
    NgivenG = Nmask - NgivenL           ; % number of gas    knowns remaining
    iMaskL  = 1 : NgivenL               ;
    iMaskG  = (NgivenL + 1) : Nmask     ;

    % Assign values for liquid-given equations
    delLL = delGiven(iMaskL)    ;
    delGL = dels    (iMaskL)    ;
    tauLL = taus    (iMaskL)    ;

    % Assign values for gas-given equations
    delLG = dels    (iMaskG)    ;
    delGG = delGiven(iMaskG)    ;
    tauGG = taus    (iMaskG)    ;
    
    % Since empty masks create row-empties, this ensures they are
    % column-empties.
    delLL = delLL(:);
    delGL = delGL(:);
    tauLL = tauLL(:);
    delLG = delLG(:);
    delGG = delGG(:);
    tauGG = tauGG(:);

    % Integer masks for vectorized Helmholtz functions:
    %       The most expensive evaluation in almost all of these thermodynamic calls is 
    %       the Helmholtz free energy functions.  As such, one the biggest optimizations
    %       lies in reducing calls to those functions by vectorizing as much as possible.
    %
    % Pack the dels and taus for input into the HFE functions
    delHelm = [delLL;delGL;delLG;delGG];
    tauHelm = [tauLL;tauLL;tauGG;tauGG];

    % Call the required HFE functions
    PhiR    = HelmholtzResidual   (delHelm,tauHelm);
    PhiR_d  = HelmholtzResidual_d (delHelm,tauHelm);
    PhiR_t  = HelmholtzResidual_t (delHelm,tauHelm);
    PhiR_dd = HelmholtzResidual_dd(delHelm,tauHelm);
    PhiR_dt = HelmholtzResidual_dt(delHelm,tauHelm);
    
    % Unpack values
    Chunks = [NgivenL,NgivenL,NgivenG,NgivenG];
    [PhiRLL   ,PhiRGL   ,PhiRLG   ,PhiRGG   ] = VectorChunk(PhiR   ,Chunks);
    [PhiRLL_d ,PhiRGL_d ,PhiRLG_d ,PhiRGG_d ] = VectorChunk(PhiR_d ,Chunks);
    [PhiRLL_t ,PhiRGL_t ,PhiRLG_t ,PhiRGG_t ] = VectorChunk(PhiR_t ,Chunks);
    [~        ,PhiRGL_dd,PhiRLG_dd,~        ] = VectorChunk(PhiR_dd,Chunks);
    [PhiRLL_dt,PhiRGL_dt,PhiRLG_dt,PhiRGG_dt] = VectorChunk(PhiR_dt,Chunks);
    

    % ============================================================================ %
    %                      Calculate Dimensionless Pressure                        %
    % ============================================================================ %
    % Pnd is the dimensionless pressure eliminated from
    % the system to shrink the solution space.  Here are the functions and
    % derivatives that are needed in the residuals below.

    fL = PhiRLL - PhiRGL + log(delLL./delGL);
    fG = PhiRLG - PhiRGG + log(delLG./delGG);
    
    R1 = [ fL .* delLL .* delGL  + delLL .* (delGL - delLL) .* (1 + delLL .* PhiRLL_d);...
           fG .* delLG .* delGG  + delLG .* (delGG - delLG) .* (1 + delLG .* PhiRLG_d)];

	R2 = [ fL .* delLL .* delGL  + delGL .* (delGL - delLL) .* (1 + delGL .* PhiRGL_d);...
           fG .* delLG .* delGG  + delGG .* (delGG - delLG) .* (1 + delGG .* PhiRGG_d)];
    
    % Helper variables for Jacobian evalutions
    fL_d = delLL .* (fL - 1 - delGL .* PhiRGL_d);
    fG_d = delGG .* (fG + 1 + delLG .* PhiRLG_d);
    fL_t = (PhiRLL_t - PhiRGL_t) .* delLL .* delGL ;
    fG_t = (PhiRLG_t - PhiRGG_t) .* delLG .* delGG ;
    %
    gLL   = 1 + delLL .* PhiRLL_d;
    gGL   = 1 + delGL .* PhiRGL_d;
    gGL_d = PhiRGL_d + delGL .* PhiRGL_dd;
    gLL_t = delLL .* PhiRLL_dt;
    gGL_t = delGL .* PhiRGL_dt;
    gLG   = 1 + delLG .* PhiRLG_d;
    gGG   = 1 + delGG .* PhiRGG_d;
    gLG_d = PhiRLG_d + delLG .* PhiRLG_dd;
    gLG_t = delLG .* PhiRLG_dt;
    gGG_t = delGG .* PhiRGG_dt;
    
    
    % Jacobian values
    R1_d = [                          delLL .* gLL + fL_d                           ;...
            (delGG - 2*delLG) .* gLG + (delGG - delLG) .* delLG .* gLG_d + fG_d    ];
        
    R1_t = [delLL .* (delGL - delLL) .* gLL_t + fL_t   ;...
            delLG .* (delGG - delLG) .* gLG_t + fG_t   ];

	R2_d = [(2*delGL - delLL) .* gGL + (delGL - delLL) .* delGL .* gGL_d + fL_d    ;...
                        -delGG .* gGG + fG_d                          ];

    R2_t = [delGL .* (delGL - delLL) .* gGL_t + fL_t   ;...
            delGG .* (delGG - delLG) .* gGG_t + fG_t   ];

    % Determinant
    DetR = R1_d .* R2_t - R1_t.*R2_d                                            ;

    % Final Newton Directions
    ddel  =  (R1 .* R2_t - R2 .* R1_t) ./ DetR ;
    dtau  =  (R2 .* R1_d - R1 .* R2_d) ./ DetR ;

    % Final outputs
    dUnknowns = [ddel,dtau,0*ddel]  ;
    RNorm     = abs(R1) + abs(R2)   ; % L_1 norm
    
end




