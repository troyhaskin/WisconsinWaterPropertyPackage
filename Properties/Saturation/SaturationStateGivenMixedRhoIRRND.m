function [Pnd,tau,delL,delG,x] = SaturationStateGivenMixedRhoIRRND(delta,iND,varargin)
    
    %   Input patterns:
    %       o SaturationStateGivenMixedRhoIRRND(delMix,iND,PhaseCheck,iMixND)
    %       o SaturationStateGivenMixedRhoIRRND(delMix,iND,PhaseCheck,iMixND,tauSat,iSatND)
    %       o SaturationStateGivenMixedRhoIRRND(delMix,iND,PhaseCheck,iMixND,tauSat,iSatND,tauGuess)
    %
    %   For this function, properties explicitly coded as <S/sat> refer to purely saturated
    %   (x = 0 or 1 only) properties on the delMix isochor (delSat = delMix).  While all
    %   tau greater than tauSat are also in saturation, tauSat is the smallest tau value
    %   at which a two-phase state can exist (while the largest is the triple-line state);
    %   therefore, the <S/sat> values allow the secant method to bounded/guarded.
    %
    
    
    % Calculate the mixture internal energy on the triple line
    [delLt,delGt] = TriplePointDensitiesR()                 ;
    [iLNDt,iGNDt] = TriplePointInternalEnergiesND()         ;
    taut          = TriplePointTau()                        ;
    xt            = QualityFromDensity(delta,delLt,delGt)   ;
    iNDt          = iLNDt + xt.*(iGNDt - iLNDt)             ;
    
    
    % Default/optional arguments
    OptionalArguments = {true,[],[],[]}             ;   % Allocate with defaults
    nOptionalIn       = length(varargin)            ;   % Number of optionals passed
    OptionalArguments(1:nOptionalIn) = varargin     ;   % Overwrite defaults
    [PhaseCheck,tauSat,iSatND,tauGuess] = OptionalArguments{:}; %#ok<ASGLU>
    
    
    
    %   Allocation
    tauHi = delta * 0   ;
    tauLo = tauHi       ;
    iHi   = tauHi       ;
    iLo   = tauHi       ;
    tau   = tauHi       ;
    Pnd   = tauHi       ;
    delL  = tauHi       ;
    delG  = tauHi       ;
    x     = tauHi       ;
    
    
    % ============================================ %
    %               Exact Phase Check              %
    % ============================================ %
    
    %   The two regions of the Vapor dome
    lessThanTriple = delta < delLt      ;
    inDMVS         = not(lessThanTriple);
    
    if any(lessThanTriple)
        
        %   non-DMVS region taus
        [~,tauLo(lessThanTriple),~,~] = SaturationStateGivenDeltaRRND(delta(lessThanTriple));
        tauHi(lessThanTriple)         = taut                                                ;
        
        %   non-DMVS region iND
        iLo(lessThanTriple) = InternalEnergyOneRND(delta(lessThanTriple),tauLo(lessThanTriple));
        iHi(lessThanTriple) = iNDt(lessThanTriple)                                             ;
    end
    
    
    if any(inDMVS)
        
        %   DMVS region taus
        deltaDMVS             = delta(inDMVS)                                                   ;
        tauBotTop             = estimateTauDMVSRegion(deltaDMVS)                                ;
        [~,tauNearTriple,~,~] = SaturationStateGivenDeltaRRND([deltaDMVS;deltaDMVS],tauBotTop)  ;
        tauLo(inDMVS)         = tauNearTriple(   1   :end/2)                                    ;
        tauHi(inDMVS)         = tauNearTriple(end/2+1: end )                                    ;
        
        %   DMVS region iND
        iLoHi       = InternalEnergyOneRND([deltaDMVS;deltaDMVS],[tauLo(inDMVS);tauHi(inDMVS)]) ;
        iLo(inDMVS) = iLoHi(   1   :end/2)                                                      ;
        iHi(inDMVS) = iLoHi(end/2+1: end )                                                      ;
        
    end
    
    %   Check upper bound
    isSaturated = (iLo - iND) > eps();
    
    %   Check lower bound
    belowDMVS   = ((iHi - iND) > eps()) & inDMVS;
    isSaturated = isSaturated & not(belowDMVS)  ;
    
    
    % If no tauGuess is given, calculate it from a secant line.
    if isempty(tauGuess)
        tauGuess = (tauHi - tauLo)./(iHi - iLo) .* (iND - iLo) + tauLo;
    end
    
    
    % ============================================ %
    %                Two Phase Solve               %
    % ============================================ %
    if any(isSaturated)
        
        %   Help variable
        iwork = iND(isSaturated)    ;

        
        %   Calculate guess internal energy
        iMix = LocalMixtureInternalEnergy(delta(isSaturated),tauGuess(isSaturated),[],[]);
        
        %   Make one more guess (this eliminates the end points to improve convergence)
        tauk     = [tauLo(isSaturated),tauGuess(isSaturated),tauHi(isSaturated)]    ;
        ik       = [iLo(isSaturated),iMix,iHi(isSaturated)]                         ;
        Rk       = ik - iwork(:,[1,1,1])                                            ;
        tauGuess = (tauk(:,1).*(Rk(:,2)<0) + tauk(:,2) + tauk(:,3).*(Rk(:,2)>=0))/2 ;
        tauk     = [tauk,tauGuess]                                                  ;
        iMix     = LocalMixtureInternalEnergy(delta(isSaturated),tauGuess,[],[])    ;
        Rk       = [ik,iMix] - iwork(:,[1,1,1,1])                                   ;
        
        %   Order the residuals
        Sign   = sign(Rk)                                   ;
        [Rk,I] = sort(abs(Rk),2)                            ;
        n      = nnz(isSaturated)                           ;
        I      = sub2ind(size(tauk),[1:n;1:n;1:n;1:n]',I)   ;
        Rk     = Rk .* Sign(I)                              ;
        tauk   = tauk(I)                                    ;
        
        %   Solve for tau
        [Pnd(isSaturated),tau(isSaturated),delL(isSaturated),delG(isSaturated)] = ...
            solveTwoPhaseSystem(tauk(:,1:3),Rk(:,1:3),delta(isSaturated),iwork);
        
        %   Solve for other variables
        %         [Pnd(isSaturated),delL(isSaturated),delG(isSaturated)] = ...
        %             SaturationStateGivenTauRRND(tau(isSaturated));
        x(isSaturated) = QualityFromDensity(delta(isSaturated),delL(isSaturated),delG(isSaturated));
    end
    
    
    
    % ============================================ %
    %                One Phase Solve               %
    % ============================================ %
    notSaturated = not(isSaturated);
    if any(notSaturated)
        
        work = delta(notSaturated);
        
        %   Get a one-phase tau guess
        tauGuess = (1-belowDMVS).*tauLo + belowDMVS.*tauHi  ;
        tauGuess = tauGuess(notSaturated)                   ;
        
        %   Call one-phase solver
        tau (notSaturated) = TemperatureOneRRND(work,iND(notSaturated),tauGuess);
        delL(notSaturated) = work                                               ;
        delG(notSaturated) = work                                               ;
        Pnd (notSaturated) = PressureOneRND(work,tau(notSaturated))             ;
        x   (notSaturated) = NaN                                                ;
    end
    
    
end


function [Pnd,tau,delL,delG] = solveTwoPhaseSystem(tauk,Rk,delta,iND)
    
    % Iteration Setup
    tolerance     = 1E-13   ;
    iterMax       = 50      ;
    notDone       = true    ;
    iter          = 0       ;
    InotConverged = 1:length(delta);
    
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
    taukp1  = tau       ;
    delLkp1 = tau       ;
    delGkp1 = tau       ;
    iNDmix  = tau       ;
    Pndkp1  = tau       ;
    
    %   Check for zeroth-guess convergence
    convergedStepWise = abs(tauk-taukm1) < tolerance;
    convergedResidual = abs(Rk)          < tolerance;
    stagnantResidual  = abs(Rk - Rkm1)   < eps()    ;
    Converged         = convergedStepWise | convergedResidual | stagnantResidual;
    NotConverged      = not(Converged);
    
    %   Update converged index array
    Ipush         = InotConverged(Converged)    ;
    InotConverged = InotConverged(NotConverged) ;
    
    
    if any(Converged)
        
        %   Update any pre-converged values
        tau(Ipush)    = tauk(Converged) ;
        
        %   Update non-tau values
        [~,Pnd(Ipush),delL(Ipush),delG(Ipush)] = ...
            LocalMixtureInternalEnergy(delta(Ipush),tau(Ipush),[],[]);
        
        % Contract any pre-converged values
        taukm2  = taukm2 (NotConverged) ;
        taukm1  = taukm1 (NotConverged) ;
        tauk    = tauk   (NotConverged) ;
        taukp1  = tauk                  ;
        Rkm2    = Rkm2   (NotConverged) ;
        Rkm1    = Rkm1   (NotConverged) ;
        Rk      = Rk     (NotConverged) ;
        iND     = iND    (NotConverged) ;
        delta   = delta  (NotConverged) ;
        iNDmix  = iNDmix (NotConverged) ;
        
    end
    
    %   Create interpolation function handles
    invQuadInterp = @(varargin) discreteHalleyUpdate(varargin{:})   ;
    secUpdate     = @(varargin) secantUpdate(varargin{:})           ;
    
    
    % Enter Iterative loop
    while notDone
        
        %   Form filter
        updateSecant           = (abs(Rk) > 1) | (abs(tauk-taukm2) < eps())   ;
        updateInverseQuadratic = not(updateSecant)                              ;
        
        
        %         (taukm2,taukm1,tauk,Rkm2,Rkm1,Rk)
        %   Compute next iterate
        taukp1(updateInverseQuadratic) = ...
            evaluateWithFilter(invQuadInterp,updateInverseQuadratic,taukm2,taukm1,tauk,Rkm2,Rkm1,Rk);
        taukp1(updateSecant) = ...
            evaluateWithFilter(secUpdate,updateSecant,taukm1,tauk,Rkm1,Rk);
        
        %   Residual calculation
        mask = (abs(Rk) > 1E-13);
        if any(mask)
            [iNDmix(mask),Pndkp1(mask),delLkp1(mask),delGkp1(mask)] = ...
                LocalMixtureInternalEnergy(delta(mask),taukp1(mask),[],[]);
        end
        mask = not(mask);
        if any(mask)
            [iNDmix(mask),Pndkp1(mask),delLkp1(mask),delGkp1(mask)] = ...
                LocalMixtureInternalEnergy(delta(mask),taukp1(mask),delLkp1(mask),delGkp1(mask));
        end
        Rkp1   = iNDmix - iND;
        
        %   Convergence check
        convergedStepWise = abs(taukp1-tauk) < tolerance                            ;
        convergedResidual = abs(Rkp1)        < tolerance                            ;
        stagnantResidual  = abs(Rkp1 - Rk)   < eps()                                ;
        Converged         = convergedStepWise | convergedResidual | stagnantResidual;
        NotConverged      = not(Converged);
        
        % Update index array
        Ipush         = InotConverged(Converged);
        InotConverged = InotConverged(NotConverged);
        
        
        % Push converged values
        tau(Ipush)  = taukp1(Converged) ;
        Pnd(Ipush)  = Pndkp1(Converged) ;
        delL(Ipush) = delLkp1(Converged);
        delG(Ipush) = delGkp1(Converged);
        
        % Contract all other vectors
        taukm2  = taukm1 (NotConverged);
        taukm1  = tauk   (NotConverged);
        tauk    = taukp1 (NotConverged);
        taukp1  = tauk                 ;
        Rkm2    = Rkm1   (NotConverged);
        Rkm1    = Rk     (NotConverged);
        Rk      = Rkp1   (NotConverged);
        iND     = iND    (NotConverged);
        delta   = delta  (NotConverged);
        delLkp1 = delLkp1(NotConverged);
        delGkp1 = delGkp1(NotConverged);
        iNDmix  = delta;
        Pndkp1  = delta;
        
        % Loop break check
        iter    = iter + 1;
        notDone = any(NotConverged) && (iter < iterMax);
        
    end
    disp(iter);
    
end

function taukp1 = secantUpdate(taukm1,tauk,Rkm1,Rk)
    %   Secant method
    dx     = Rk .* (tauk - taukm1) ./ (Rk - Rkm1)   ;
    taukp1 = tauk - dx                              ;
    
    %   Guard
    mask = (taukp1 > TriplePointTau()) | (taukp1 < 1)   ;
    alpha       = 0.5                                   ;
    while any(mask)
        taukp1(mask) = tauk(mask) - alpha*dx(mask)                  ;
        alpha        = 0.5*alpha                                    ;
        mask         = (taukp1 > TriplePointTau()) | (taukp1 < 1)   ;
    end
    
    if any(isnan(taukp1))
        g = [ ];
    end
    
end


function taukp1 = discreteHalleyUpdate(taukm2,taukm1,tauk,Rkm2,Rkm1,Rk)
    %   Secant method
    idtaukkm1   = 1./(tauk   - taukm1);
    idtaukkm2   = 1./(tauk   - taukm2);
    idtaukm1km2 = 1./(taukm1 - taukm2);
    
    DRk  = (idtaukkm1 + idtaukkm2).*Rk  - (idtaukkm1 + idtaukm1km2).*Rkm1 + (idtaukm1km2 - idtaukkm2).*Rkm2;
    DDRk = 2*((idtaukkm1 .* idtaukkm2).*Rk - (idtaukkm1 .* idtaukm1km2).*Rkm1 + (idtaukkm2 .* idtaukm1km2).*Rkm2);
    
    dx     = 2*Rk.*DRk./ (2*DRk.^2 - Rk.* DDRk) ;
    taukp1 = tauk - dx                          ;
    
    mask = (taukp1 > TriplePointTau()) | (taukp1 < 1);
    alpha       = 0.5                                       ;
    while any(mask)
        taukp1(mask) = tauk(mask) - alpha*dx(mask)                  ;
        alpha        = 0.5*alpha                                    ;
        mask         = (taukp1 > TriplePointTau()) | (taukp1 < 1)   ;
    end
    
    if any(isnan(taukp1))
        g = [ ];
    end
    
end




function [iNDmix,Pnd,delL,delG] = LocalMixtureInternalEnergy(del,tau,delL,delG)
    
    % Get the saturation state
    [Pnd,delL,delG] = SaturationStateGivenTauRRND(tau,delL,delG);
    
    % Calculate the mixture internal energy
    x     = QualityFromDensity(del,delL,delG);
    iLGND = InternalEnergyOneRND([delL;delG],[tau;tau]);
    iLND  = iLGND(1:end/2)      ;
    iGND  = iLGND(end/2+1:end)  ;
    iNDmix = iLND + x.*(iGND - iLND);
    
end





