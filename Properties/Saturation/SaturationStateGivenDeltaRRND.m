function [Pnd,tau,delL,delG] = SaturationStateGivenDeltaRRND(delta,tau0)
    
    % ========================================================= %
    %                        Set-Up                             %
    % ========================================================= %
    
    % Look for given guess tau
    if (nargin < 2) || isempty(tau0)
        useDeltaEstimation = true;
    else
        useDeltaEstimation = false;
    end

    % Given Masks
    givenL = delta >  1 ;
    givenG = delta <  1 ;
    givenC = delta == 1 ;

    % Ordering Masks: since the solution algorithm relies on the ordering of dels to be
    % [delL;delG], these index arrays allow insertion of the solved properties into the 
    % original ordering.
    nGiven    = length(delta);
    Ioriginal = 1 : nGiven;
    Iordered  = [Ioriginal(givenL),Ioriginal(givenG),Ioriginal(givenC)];

    % Given densities
    delLgiven = delta(givenL)                   ;
    delGgiven = delta(givenG)                   ;
    delCgiven = delta(givenC)                   ;
    delGiven  = [delLgiven;delGgiven;delCgiven] ; % destroys original ordering; use index arrays above

    % Number of gas densities given
    nGivenL = nnz(givenL);
    nGivenG = nnz(givenG);
    nGivenC = nnz(givenC);


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
    [delLnearC , delGnearC] = NearCriticalStabilityDensityLimitsR() ;
    [delLt,~]               = TriplePointDensitiesR()               ;
    [delGmin,delLmax]       = saturableDeltas()                     ;

    %   Density is within the bounds of the vapor dome.
    canSaturate    = [delLgiven <= delLmax ; delGgiven >= delGmin ; true(nGivenC,1)];
    cannotSaturate = not(canSaturate)                                                       ;
    
    %   Near critical densities are avoided due to stability issues (near singular Jacobian)
    shouldIterate = [delLgiven >= delLnearC ; delGgiven <= delGnearC ; false(nGivenC,1)];   

    %   Density is within the DMVS region
    nearTriple = [(delLgiven >= delLt) & (delLgiven <= delLmax) ; false(nGivenG+nGivenC,1)];

    % Logical index for values that will be iterated upon
    willGuess   = (not(nearTriple) | not(useDeltaEstimation)) & canSaturate ;
    willIterate = willGuess & shouldIterate                                 ;



    % ================================================================================== %
    %                                    Guess Values                                    %
    % ================================================================================== %
    
    %   Allocation
    delLguess = delta;
    delGguess = delta;
    tauGuess  = delta;


    % Get guess values for all given density values not near the triple line
    if useDeltaEstimation
        
        %   No tau0 given
        if any(willGuess)
            [delLguess,delGguess,tauGuess] = ...
                AssignWithFilter(...
                    @(wg) EstimateDelLDelGTauFromDel(delGiven(wg)),...
                    willGuess,delLguess,delGguess,tauGuess);
        end

    else

         if any(canSaturate)
            %   tau0 given
            tauGuess = tau0;
            delLguess(givenG & canSaturate) = EstimateDelLFromTau(tauGuess(givenG & canSaturate));
            delGguess(givenL & canSaturate) = EstimateDelGFromTau(tauGuess(givenL & canSaturate));
        end

    end


    % ================================================================================== %
    %                                 Newton Solution                                    %
    % ================================================================================== %
    
    if any(willIterate)

        % Iteration parameters
        Tolerance   = DefaultAbsoluteIterationTolerance();
        IterMax     = DefaultMaximumIterationCount();

        % Fill guess matrix
        xSol = [[delGguess(1:nGivenL);delLguess(nGivenL+1:end)] , tauGuess];
        xSol = xSol(willIterate,:);
        nLIterate = nnz(delGiven(willIterate) > 1);
        
        %   Perform one Newton calculation to initialize Broyden Jacobian
        [dx,~,S.r1,S.r2,S.iJ11,S.iJ12,S.iJ21,S.iJ22] = ...
            newton(xSol,1:nnz(willIterate),delGiven(willIterate),nLIterate);
        S.del = xSol(:,1);
        S.tau = xSol(:,2);

        % Newton solution
        %   Solve
        xSol = NewtonUpdater(@(x,mask) broydenClose(x,mask,delGiven(willIterate),nLIterate),xSol-dx,1E-13,IterMax);
        xSol = NewtonUpdater(@(x,mask)       newton(x,mask,delGiven(willIterate),nLIterate),xSol,Tolerance,IterMax);
        
        
    else
        xSol = delGiven(:,[1,1])  ;
    end
    
    % Allocations
    delG = delGguess ;
    delL = delLguess ;
    tau  = tauGuess  ;
    
    % Insert Newton solutions
    Iwill                          = Ioriginal(willIterate)         ;
    delG(Iwill(1:nLIterate))       = xSol(1:nLIterate,1)        ;
    delL(Iwill((nLIterate+1):end)) = xSol((nLIterate+1):end,1)  ;
    tau (Iwill)                    = xSol(   :   ,2)            ;

    %   Critical quantities
    delG(givenC) = 1;
    delL(givenC) = 1;
    tau (givenC) = 1;
    
    % Undefinable quantities
    delG(cannotSaturate) = 0  ;
    delL(cannotSaturate) = 0  ;
    tau (cannotSaturate) = 0  ;

    % Put back into passed-in ordering
    delG(Iordered) = delG;
    delL(Iordered) = delL;
    tau (Iordered) = tau;

    % Get Pressures
    Pnd = PressureOneRND(delG,tau);








    function [dx,rNorm] = broydenClose(x,mask,delGiven,NliquidGiven)
        [dx,rNorm,S] = broyden(x,mask,delGiven,NliquidGiven,S);
    end

end



function [dx,rNorm,r1,r2,iJ11,iJ12,iJ21,iJ22] = newton(x,mask,delGiven,nLiquidGiven)
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
    dels     = x(:,1)    ;
    taus     = x(:,2)    ;
    
    %   Indexes/masks for determining liquid vs. gas given
    Nmask   = length(mask)              ;
    NgivenL = nnz(mask <= nLiquidGiven) ; % number of liquid knowns remaining
    NgivenG = Nmask - NgivenL           ; % number of gas    knowns remaining
    iMaskL  =       1     : NgivenL     ;
    iMaskG  = (NgivenL+1) :  Nmask      ;

    % Assign values for liquid-given equations
    delLL = delGiven(mask(iMaskL))  ;
    delGL = dels    (iMaskL)        ;
    tauLL = taus    (iMaskL)        ;

    % Assign values for gas-given equations
    delLG = dels    (iMaskG)        ;
    delGG = delGiven(mask(iMaskG))  ;
    tauGG = taus    (iMaskG)        ;
    
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
    [PhiR,PhiR_d,PhiR_t,PhiR_dt,PhiR_dd] = HelmholtzResidualCombo__d_t_dt_dd(delHelm,tauHelm);
    
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
    
    r1 = [ fL .* delLL .* delGL  + delLL .* (delGL - delLL) .* (1 + delLL .* PhiRLL_d);...
           fG .* delLG .* delGG  + delLG .* (delGG - delLG) .* (1 + delLG .* PhiRLG_d)];

	r2 = [ fL .* delLL .* delGL  + delGL .* (delGL - delLL) .* (1 + delGL .* PhiRGL_d);...
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
    r1_d = [                          delLL .* gLL + fL_d                           ;...
            (delGG - 2*delLG) .* gLG + (delGG - delLG) .* delLG .* gLG_d + fG_d    ];
        
    r1_t = [delLL .* (delGL - delLL) .* gLL_t + fL_t   ;...
            delLG .* (delGG - delLG) .* gLG_t + fG_t   ];

	r2_d = [(2*delGL - delLL) .* gGL + (delGL - delLL) .* delGL .* gGL_d + fL_d    ;...
                        -delGG .* gGG + fG_d                          ];

    r2_t = [delGL .* (delGL - delLL) .* gGL_t + fL_t   ;...
            delGG .* (delGG - delLG) .* gGG_t + fG_t   ];

    % Determinant
    detJ = r1_d .* r2_t - r1_t.*r2_d;

    
    % Inverse Jacobian elements
    iJ11 =  r2_t ./ detJ;
    iJ12 = -r1_t ./ detJ;
    iJ21 = -r2_d ./ detJ;
    iJ22 =  r1_d ./ detJ;

    % Final Newton Directions
    ddel  =  iJ11 .* r1 + iJ12 .* r2;
    dtau  =  iJ21 .* r1 + iJ22 .* r2;

    % Final outputs
    dx    = [ddel,dtau]             ;
    rNorm = abs(real(r1)) + abs(real(r2)) + imag(r1) + imag(r2); % L_1 norm
    
end



function [dx,rNorm,S] = broyden(x,mask,delGiven,nLgiven,S)
    
    % Pull unknowns
    del = x(:,1)    ;
    tau = x(:,2)    ;
    
    %   Indexes/masks for determining liquid vs. gas given
    Nmask   = length(mask)          ;
    NgivenL = nnz(mask <= nLgiven)  ; % number of liquid knowns remaining
    NgivenG = Nmask - NgivenL       ; % number of gas    knowns remaining
    iMaskL  =       1     : NgivenL ;
    iMaskG  = (NgivenL+1) :  Nmask  ;

    % Assign values for liquid-given equations
    delLL = delGiven(mask(iMaskL)) ;
    delGL = del    (iMaskL)        ;
    tauLL = tau    (iMaskL)        ;

    % Assign values for gas-given equations
    delLG = del    (iMaskG)        ;
    delGG = delGiven(mask(iMaskG))  ;
    tauGG = tau    (iMaskG)        ;
    
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
    [PhiR,PhiR_d] = HelmholtzResidualCombo__d(delHelm,tauHelm);
    
    % Unpack values
    Chunks = [NgivenL,NgivenL,NgivenG,NgivenG];
    [PhiRLL   ,PhiRGL   ,PhiRLG   ,PhiRGG   ] = VectorChunk(PhiR   ,Chunks);
    [PhiRLL_d ,PhiRGL_d ,PhiRLG_d ,PhiRGG_d ] = VectorChunk(PhiR_d ,Chunks);
    

    % ============================================================================ %
    %                      Calculate Dimensionless Pressure                        %
    % ============================================================================ %
    % Pnd is the dimensionless pressure eliminated from
    % the system to shrink the solution space.  Here are the functions and
    % derivatives that are needed in the residuals below.

    fL = PhiRLL - PhiRGL + log(delLL./delGL);
    fG = PhiRLG - PhiRGG + log(delLG./delGG);
    
    r1 = [ fL .* delLL .* delGL  + delLL .* (delGL - delLL) .* (1 + delLL .* PhiRLL_d);...
           fG .* delLG .* delGG  + delLG .* (delGG - delLG) .* (1 + delLG .* PhiRLG_d)];

	r2 = [ fL .* delLL .* delGL  + delGL .* (delGL - delLL) .* (1 + delGL .* PhiRGL_d);...
           fG .* delLG .* delGG  + delGG .* (delGG - delLG) .* (1 + delGG .* PhiRGG_d)];

       
     % Form deltas
    dx1  = del - S.del(mask)    ;
    dx2  = tau - S.tau(mask)    ;
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
    ddel = S.iJ11(mask) .* r1 + S.iJ12(mask) .* r2;
    dtau = S.iJ21(mask) .* r1 + S.iJ22(mask) .* r2;
    
    % Pack updates and calculate norm for the Newton updater
    dx    = [ddel,dtau]         ;
    rNorm = abs(r1) + abs(r2)   ;
    
    
    %   Update struct elements
    S.del(mask) = del  ;
    S.tau(mask) = tau  ;
    S.r1(mask)  = r1   ;
    S.r2(mask)  = r2   ;
end

