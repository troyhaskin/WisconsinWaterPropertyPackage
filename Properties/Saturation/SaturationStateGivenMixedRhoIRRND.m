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
    isSaturated = iND < iLo ;
    
    %   Check lower bound
    belowDMVS   = (iND < iHi) & inDMVS          ;
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
%         dwork = delta(isSaturated)  ;
        
        %   Define the three interpolation points used for two-phase solution
        tauk   = [tauLo(isSaturated),tauGuess(isSaturated),tauHi(isSaturated)];
        iMix   = LocalMixtureInternalEnergy(delta(isSaturated),tauk(:,2),[],[])     ;
        Rk     = [iLo(isSaturated),iMix,iHi(isSaturated)] - iwork(:,[1,1,1]);
        [Rk,I] = sort(Rk,2);
        n      = nnz(isSaturated);
        tauk   = tauk(sub2ind(size(tauk),[1:n;1:n;1:n]',I));
        
        %   Solve for tau
        [Pnd(isSaturated),tau(isSaturated),delL(isSaturated),delG(isSaturated)] = ...
            solveTwoPhaseSystem(tauk(:,1),tauk(:,2),tauk(:,3),Rk(:,1),Rk(:,2),Rk(:,3),...
            delta(isSaturated),iwork);
        
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


function [Pnd,tau,delL,delG] = solveTwoPhaseSystem(taukm2,taukm1,tauk,Rkm2,Rkm1,Rk,delta,iND)
    
    % Iteration Setup
    tolerance     = 1E-13   ;
    iterMax       = 50      ;
    notDone       = true    ;
    iter          = 1       ;
    InotConverged = 1:length(delta);
    
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
    Converged    = abs(Rk) < tolerance;
    NotConverged = not(Converged);
    
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
    invQuadInterp = @(varargin) discreteHalleyUpdate(varargin{:})  ;
    secUpdate     = @(varargin) secantUpdate(varargin{:})                   ;

    
    % Enter Iterative loop
    while notDone
        
        %   Form filter
        updateSecant           = abs(Rk) > 1    ;
        updateInverseQuadratic = not(updateSecant)  ;
        
        
%         (taukm2,taukm1,tauk,Rkm2,Rkm1,Rk)
        %   Compute next iterate
        taukp1(updateInverseQuadratic) = ...
            evaluateWithFilter(invQuadInterp,updateInverseQuadratic,taukm2,taukm1,tauk,Rkm2,Rkm1,Rk);
        taukp1(updateSecant) = ...
            evaluateWithFilter(secUpdate,updateSecant,taukm1,tauk,Rkm1,Rk);

        %   Residual calculation
        mask = abs(Rk) > 1E-13;
        [iNDmix(mask),Pndkp1(mask),delLkp1(mask),delGkp1(mask)] = ...
                LocalMixtureInternalEnergy(delta(mask),taukp1(mask),[],[]);
        mask = not(mask);
        [iNDmix(mask),Pndkp1(mask),delLkp1(mask),delGkp1(mask)] = ...
                LocalMixtureInternalEnergy(delta(mask),taukp1(mask),delLkp1(mask),delGkp1(mask));
        Rkp1   = iNDmix - iND;
        
        %   Convergence check
        Converged    = (abs(taukp1-tauk) < tolerance) | (abs(Rkp1) < tolerance) ;
        NotConverged = not(Converged);
        
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
    taukp1 = tauk - Rk .* (tauk - taukm1) ./ (Rk - Rkm1);
    
    mask = (taukp1 > TriplePointTau()) | (taukp1 < 1);
    alpha       = 0.5                                       ;
    while any(mask)
        dx           = Rk(mask).*(tauk(mask)-taukm1(mask))./(Rk(mask)-Rkm1(mask))   ;
        taukp1(mask) = tauk(mask) - alpha*dx                                        ;
        alpha        = 0.5*alpha                                                    ;
        mask         = (taukp1 > TriplePointTau()) | (taukp1 < 1)                   ;
    end

end

function taukp1 = discreteHalleyUpdate(taukm2,taukm1,tauk,Rkm2,Rkm1,Rk)
    %   Secant method
    idtaukkm1   = 1./(tauk   - taukm1 + eps());
    idtaukkm2   = 1./(tauk   - taukm2 + eps());
    idtaukm1km2 = 1./(taukm1 - taukm2 + eps());
    
    DRk  = (idtaukkm1 + idtaukkm2).*Rk  - (idtaukkm1 + idtaukm1km2).*Rkm1 + (idtaukm1km2 - idtaukkm2).*Rkm2;
    DDRk = 2*((idtaukkm1 .* idtaukkm2).*Rk - (idtaukkm1 .* idtaukm1km2).*Rk + (idtaukkm2 .* idtaukm1km2).*Rk);
    
    taukp1 = tauk - 2*Rk.*DRk./ (2*DRk.^2 - Rk.* DDRk);
    
    mask = (taukp1 > TriplePointTau()) | (taukp1 < 1);
    alpha       = 0.5                                       ;
    while any(mask)
        dx           = 2*Rk(mask).*DRk(mask)./ (2*DRk(mask).^2 - Rk(mask).* DDRk(mask)) ;
        taukp1(mask) = tauk(mask) - alpha*dx                                            ;
        alpha        = 0.5*alpha                                                        ;
        mask         = (taukp1 > TriplePointTau()) | (taukp1 < 1)                       ;
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





