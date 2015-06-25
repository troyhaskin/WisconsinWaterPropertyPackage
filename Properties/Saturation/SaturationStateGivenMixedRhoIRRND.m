function state = SaturationStateGivenMixedRhoIRRND(delta,iND,tauGuess)
    
    %   Input patterns:
    %       o SaturationStateGivenMixedRhoIRRND(delMix,iND,tauGuess)
    %
    %   For this function, properties explicitly coded as <S/sat> refer to purely saturated
    %   (x = 0 or 1 only) properties on the delMix isochor (delSat = delMix).  While all
    %   tau greater than tauSat are also in saturation, tauSat is the smallest tau value
    %   at which a two-phase state can exist (while the largest is the triple-line state);
    %   therefore, the <S/sat> values allow the secant method to bounded/guarded.
    %
    
    % Default/optional arguments
    if (nargin < 3)
        tauGuess = [];
    end
    
    
    %   Calculate the mixture internal energy on the triple line
    [delLt,delGt] = TriplePointDeltas()                     ;
    [iLNDt,iGNDt] = TriplePointInternalEnergiesND()         ;
    taut          = TriplePointTau()                        ;
    xt            = QualityFromDensity(delta,delLt,delGt)   ;
    iNDt          = iLNDt + xt.*(iGNDt - iLNDt)             ;



    %   Allocation
    N               = numel(delta)  ;
    tauHi(N,1)      = 0             ;
    tauLo(N,1)      = 0             ;
    iHi(N,1)        = 0             ;
    iLo(N,1)        = 0             ;
    state.tau(N,1)  = 0             ;
    state.Pnd(N,1)  = 0             ;
    state.delL(N,1) = 0             ;
    state.delG(N,1) = 0             ;
    state.x(N,1)    = 0             ;

    
    
    % ============================================ %
    %               Exact Phase Check              %
    % ============================================ %

    %   Catch any points on the triple line
    isTripleTau      = abs(iNDt - iND) < eps()                                                  ;
    isCriticalTau    = abs(delta - 1)  < eps() & abs(iND - DimensioningInternalEnergy()) < eps();
    isNotTripleTau   = not(isTripleTau)                                                         ;
    isNotCriticalTau = not(isCriticalTau)                                                       ;


    %   The two regions of the Vapor dome
    lessThanTriple = delta < delLt      ;
    inDMVS         = not(lessThanTriple);
    
    
    %   Handle triple and critical states
    iLo(isTripleTau)         = iNDt(isTripleTau)            ;
    iHi(isTripleTau)         = iNDt(isTripleTau)            ;
    state.tau(isTripleTau)   = taut                         ;
    iLo(isCriticalTau)       = DimensioningInternalEnergy() ;
    iHi(isCriticalTau)       = DimensioningInternalEnergy() ;
    state.tau(isCriticalTau) = 1                            ;

    
    %   Non-DMVS region saturation calculation
    calculate = lessThanTriple & isNotTripleTau & isNotCriticalTau;
    if any(calculate)
        
        %   non-DMVS region taus
        [~,tauLo(calculate),~,~] = SaturationStateGivenDeltaRRND(delta(calculate))  ;
        tauHi(calculate)         = taut                                             ;
        
        %   non-DMVS region iND
        iLo(calculate) = InternalEnergyOneRND(delta(calculate),tauLo(calculate));
        iHi(calculate) = iNDt(calculate)                                        ;
    end
    

    %   DMVS region saturation calculation
    calculate = inDMVS & isNotTripleTau & isNotCriticalTau;
    if any(calculate)
        
        %   DMVS region taus
        deltaDMVS             = delta(calculate)                                                ;
        tauBotTop             = estimateTauDMVSRegion(deltaDMVS)                                ;
        [~,tauNearTriple,~,~] = SaturationStateGivenDeltaRRND([deltaDMVS;deltaDMVS],tauBotTop)  ;
        tauLo(calculate)      = tauNearTriple(   1   :end/2)                                    ;
        tauHi(calculate)      = tauNearTriple(end/2+1: end )                                    ;
        
        %   DMVS region iND
        iLoHi = InternalEnergyOneRND([deltaDMVS;deltaDMVS],[tauLo(calculate);tauHi(calculate)]) ;
        iLo(calculate) = iLoHi(   1   :end/2)                                                   ;
        iHi(calculate) = iLoHi(end/2+1: end )                                                   ;
        
    end


    %   Check upper bound
    isSaturated = (iLo - iND) > eps();
    
    %   Check lower bound
    belowDMVS   = ((iHi - iND) > eps()) & inDMVS;
    isSaturated = isSaturated & not(belowDMVS)  ;
    
    %   Save exact phase predicate in output
    state.isTwoPhi = isSaturated;
    
    % If no tauGuess is given, calculate it from a secant.
    if isempty(tauGuess)
        tauGuess = (tauHi - tauLo)./(iHi - iLo) .* (iND - iLo) + tauLo;
    end
    
    
    % ============================================ %
    %                Two Phase Solve               %
    % ============================================ %
    calculate = isSaturated & isNotTripleTau & isNotCriticalTau;
    if any(calculate)
        
        %   Help variable
        iwork = iND(calculate);
        
        %   Calculate guess internal energy
        iMix = LocalMixtureInternalEnergy(delta(calculate),tauGuess(calculate));
        
        %   Make one more guess (this eliminates the end points to improve convergence)
        tauk     = [tauLo(calculate),tauGuess(calculate),tauHi(calculate)]          ;
        ik       = [iLo(calculate),iMix,iHi(calculate)]                             ;
        Rk       = ik - iwork(:,[1,1,1])                                            ;
        tauGuess = (tauk(:,1).*(Rk(:,2)<0) + tauk(:,2) + tauk(:,3).*(Rk(:,2)>=0))/2 ;
        tauk     = [tauk,tauGuess]                                                  ;
        iMix     = LocalMixtureInternalEnergy(delta(calculate),tauGuess)            ;
        Rk       = [ik,iMix] - iwork(:,[1,1,1,1])                                   ;
        
        %   Order the residuals
        Sign   = sign(Rk)                               ;
        [Rk,I] = sort(abs(Rk),2)                        ;
        n      = (1:nnz(calculate))'                    ;
        I      = sub2ind(size(tauk),n(:,[1,1,1,1]),I)   ;
        Rk     = Rk .* Sign(I)                          ;
        tauk   = tauk(I)                                ;
        
        %   Solve for tau
        [state.Pnd(calculate),state.tau(calculate),state.delL(calculate),state.delG(calculate)] = ...
            solveTwoPhaseSystem(tauk(:,1:3),Rk(:,1:3),delta(calculate),iwork);
        
        %   Quality
        state.x(calculate) = QualityFromDensity(delta(calculate),state.delL(calculate),state.delG(calculate));
    end
    
    
    
    % ============================================ %
    %                One Phase Solve               %
    % ============================================ %
    calculate = not(isSaturated) & isNotTripleTau & isNotCriticalTau;
    if any(calculate)
        
        work = delta(calculate);
        
        %   Get a one-phase tau guess
        tauGuess = (1-belowDMVS).*tauLo + belowDMVS.*tauHi  ;
        tauGuess = tauGuess(calculate)                      ;
        
        %   Call one-phase solver
        state.tau (calculate) = TemperatureOneRRND(work,iND(calculate),tauGuess)  ;
        state.delL(calculate) = work                                              ;
        state.delG(calculate) = work                                              ;
        state.Pnd (calculate) = PressureOneRND(work,state.tau(calculate))         ;
        state.x   (calculate) = NaN                                               ;
    end
    
    
end


function [Pnd,tau,delL,delG] = solveTwoPhaseSystem(tauk,Rk,delta,iND)
    
    % Iteration Setup
    tolerance     = 1E-13           ;   %   Confidently achievable tolerance for double precision 
    iterMax       = 50              ;
    iter          = 0               ;
    InotConverged = 1:length(delta) ;
    
    %   Extract columns to variables
    taukm2 = tauk(:,3)  ;
    taukm1 = tauk(:,2)  ;
    tauk   = tauk(:,1)  ;
    Rkm2   = Rk(:,3)    ;
    Rkm1   = Rk(:,2)    ;
    Rk     = Rk(:,1)    ;
    

    %   Allocation
    tau     = tauk * 0  ;
    taukp1  = tau       ; %#ok<NASGU>
    Pnd     = tau       ;
    delL    = tau       ;
    delG    = tau       ;
    
    %   Check for zeroth-guess convergence
    convergedStepWise = abs(tauk-taukm1) < tolerance                            ;
    convergedResidual = abs(Rk)          < tolerance                            ;
    stagnantResidual  = abs(Rk - Rkm1)   < eps()                                ;
    converged         = convergedStepWise | convergedResidual | stagnantResidual;
    notConverged      = not(converged)                                          ;
    notDone           = any(notConverged)                                       ;
    
    %   Update converged index array
    Ipush         = InotConverged(converged)    ;
    InotConverged = InotConverged(notConverged) ;
    
    
    if any(converged)
        
        %   Update any pre-converged values
        tau(Ipush)    = tauk(converged) ;
        
        %   Update non-tau values
        [~,Pnd(Ipush),delL(Ipush),delG(Ipush)] = ...
            LocalMixtureInternalEnergy(delta(Ipush),tau(Ipush));
        
        % Contract any pre-converged values
        taukm2  = taukm2 (notConverged) ;
        taukm1  = taukm1 (notConverged) ;
        tauk    = tauk   (notConverged) ;
        Rkm2    = Rkm2   (notConverged) ;
        Rkm1    = Rkm1   (notConverged) ;
        Rk      = Rk     (notConverged) ;
        iND     = iND    (notConverged) ;
        delta   = delta  (notConverged) ;
        
    end
    
    
    % Enter Iterative loop
    while notDone

        %   Compute next iterate
        taukp1 = discreteNewtonUpdate(taukm2,taukm1,tauk,Rkm2,Rkm1,Rk);
       
        %   Residual calculation
        [iNDmix,Pndkp1,delLkp1,delGkp1] = LocalMixtureInternalEnergy(delta,taukp1);
        Rkp1   = iNDmix - iND;
        
        %   Convergence check
        convergedStepWise = abs(taukp1-tauk) < tolerance                            ;
        convergedResidual = abs(Rkp1)        < tolerance                            ;
        stagnantResidual  = abs(Rkp1 - Rk)   < eps()                                ;
        converged         = convergedStepWise | convergedResidual | stagnantResidual;
        notConverged      = not(converged);
        
        % Update index array
        Ipush         = InotConverged(converged);
        InotConverged = InotConverged(notConverged);
        
        % Push converged values
        tau(Ipush)  = taukp1(converged) ;
        Pnd(Ipush)  = Pndkp1(converged) ;
        delL(Ipush) = delLkp1(converged);
        delG(Ipush) = delGkp1(converged);
        
        % Contract all other vectors
        taukm2  = taukm1 (notConverged);
        taukm1  = tauk   (notConverged);
        tauk    = taukp1 (notConverged);
        Rkm2    = Rkm1   (notConverged);
        Rkm1    = Rk     (notConverged);
        Rk      = Rkp1   (notConverged);
        iND     = iND    (notConverged);
        delta   = delta  (notConverged);
        
        % Loop break check
        iter    = iter + 1;
        notDone = any(notConverged) && (iter < iterMax);
        
    end
    
end

function taukp1 = discreteNewtonUpdate(taukm2,taukm1,tauk,Rkm2,Rkm1,Rk)
    %   Secant method
    dx1 = tauk - taukm1;
    dx2 = tauk - taukm2;
    DRk = (1./dx1+1./dx2).*Rk + dx2./(dx1.^2-dx1.*dx2).*Rkm1+dx1./(dx2.*(dx2-dx1)).*Rkm2;
    
    dx     = Rk ./DRk               ;
    taukp1 = guardedStep(tauk,dx)   ;
    
end

function taukp1 = guardedStep(tau,dx)
    
    taukp1 = tau - dx                                   ;
    mask   = (taukp1 > TriplePointTau()) | (taukp1 < 1) ;
    alpha  = 1                                          ;

    while any(mask)
        taukp1(mask) = tau(mask) - alpha*dx(mask)                   ;
        alpha        = 0.5*alpha                                    ;
        mask         = (taukp1 > TriplePointTau()) | (taukp1 < 1)   ;
    end

end


function [iNDmix,Pnd,delL,delG] = LocalMixtureInternalEnergy(del,tau)
    
    % Get the saturation state
    [Pnd,delL,delG] = SaturationStateGivenTauRRND(tau);
    
    % Calculate the mixture internal energy
    x     = QualityFromDensity(del,delL,delG);
    iLGND = InternalEnergyOneRND([delL;delG],[tau;tau]);
    iLND  = iLGND(1:end/2)      ;
    iGND  = iLGND(end/2+1:end)  ;
    iNDmix = iLND + x.*(iGND - iLND);
    
end

