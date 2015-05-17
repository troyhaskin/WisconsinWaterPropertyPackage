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
        work = iND(isSaturated);
        
        %   Define the three interpolation points used for two-phase solution
        tauk   = [tauLo(isSaturated),tauGuess(isSaturated),tauHi(isSaturated)];
        iMix   = LocalMixtureInternalEnergy(delta(isSaturated),tauk(:,2))     ;
        Rk     = [iLo(isSaturated),iMix,iHi(isSaturated)] - work(:,[1,1,1]);
        [Rk,I] = sort(Rk,2);
        n      = nnz(isSaturated);
        tauk   = tauk(sub2ind(size(tauk),[1:n;1:n;1:n]',I));
        
        %   Solve for tau
        tau(isSaturated) = solveTwoPhaseSystem(...
            tauk(:,1),tauk(:,2),tauk(:,3),Rk(:,1),Rk(:,2),Rk(:,3),delta(isSaturated),work);

        %   Solve for other variables
        [Pnd(isSaturated),delL(isSaturated),delG(isSaturated)] = ...
            SaturationStateGivenTauRRND(tau(isSaturated));
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


function tau = solveTwoPhaseSystem(taukm2,taukm1,tauk,Rkm2,Rkm1,Rk,delta,iND)

    % Iteration Setup
    tolerance     = 1E-12 ;
    iterMax       = 1E2   ;
    notDone       = true  ;
    iter          = 0     ;
    InotConverged = 1:length(delta);
    
    % Check for zeroth-guess convergence
    Converged    = abs(Rk) < tolerance;
    NotConverged = not(Converged);
    
    % Update converged index array
    Ipush         = InotConverged(Converged)    ;
    InotConverged = InotConverged(NotConverged) ;
    tau(Ipush)    = tauk(Converged)             ;

    % Contract all other vectors
    taukm2  = taukm2(NotConverged);
    taukm1  = taukm1(NotConverged);
    tauk    = tauk  (NotConverged);
    taukp1  = tauk                ;
    Rkm2    = Rkm2  (NotConverged);
    Rkm1    = Rkm1  (NotConverged);
    Rk      = Rk    (NotConverged);
    iND     = iND   (NotConverged);
    delta   = delta (NotConverged);
    
    %   Create function handles
    quadInterp    = @(varargin) quadraticInterpolation(varargin{:})         ;
    invQuadInterp = @(varargin) inverseQuadraticInterpolation(varargin{:})  ;
    
    % Enter Iterative loop
    while notDone
    
        %   Form filter
        interpQuadratic        = abs(Rk) > 1;
        interpInverseQuadratic = not(interpQuadratic);
        
        %   Compute next iterate
        taukp1(interpQuadratic) = ...
            evaluateWithFilter(quadInterp,interpQuadratic,taukm2,taukm1,tauk,Rkm2,Rkm1,Rk);
        taukp1(interpInverseQuadratic) = ...
            evaluateWithFilter(invQuadInterp,interpInverseQuadratic,taukm2,taukm1,tauk,Rkm2,Rkm1,Rk);

        %   Residual calculation
        Rkp1   = LocalMixtureInternalEnergy(delta,taukp1) - iND;

        %   Convergence check
        Converged    = (abs(taukp1-tauk) < tolerance) | (abs(Rkp1) < tolerance);
        NotConverged = not(Converged);

        % Update index array
        Ipush         = InotConverged(Converged);
        InotConverged = InotConverged(NotConverged);

        % Push converged values
        tau(Ipush) = taukp1(Converged);

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
        
        
%     %   Interpolation
%         taukp1 = Rkm1.*Rk  ./((Rkm2-Rkm1).*(Rkm2-Rk)) .* taukm2 + ...
%                  Rkm2.*Rk  ./((Rkm1-Rkm2).*(Rkm1-Rk)) .* taukm1 + ...
%                  Rkm2.*Rkm1./((Rkm2-Rk  ).*(Rkm1-Rk)) .* tauk   ;
        
        
        
        % Loop break check
        iter    = iter + 1;
        notDone = any(NotConverged) && (iter < iterMax);

    end
disp(iter)
end

function taukp1 = quadraticInterpolation(taukm2,taukm1,tauk,Rkm2,Rkm1,Rk)
        dtaukkm1   = tauk   - taukm1                ;
        dtaukkm2   = tauk   - taukm2                ;
        dtaukm1km2 = taukm1 - taukm2                ;
        fkkm1      = (Rk   - Rkm1)./dtaukkm1        ;
        fkkm2      = (Rk   - Rkm2)./dtaukkm2        ;
        fkm1km2    = (Rkm1 - Rkm2)./dtaukm1km2      ;
        fkkm1km2   = Rk  ./(dtaukkm1.*dtaukkm2  )   -...
                     Rkm1./(dtaukkm1.*dtaukm1km2)   +...
                     Rkm2./(dtaukkm2.*dtaukm1km2)   ;
        w          = fkkm1 + fkkm2 + fkm1km2        ;
        taukp1     = tauk - 2*Rk ./(w + sign(w).*sqrt(w.^2 - 4* Rk.*fkkm1km2));
end


function taukp1 = inverseQuadraticInterpolation(taukm2,taukm1,tauk,Rkm2,Rkm1,Rk)
    %   Inverse quadratic interpolation afterward
    taukp1 = Rkm1.*Rk  ./((Rkm2-Rkm1).*(Rkm2-Rk)) .* taukm2 + ...
             Rkm2.*Rk  ./((Rkm1-Rkm2).*(Rkm1-Rk)) .* taukm1 + ...
             Rkm2.*Rkm1./((Rkm2-Rk  ).*(Rkm1-Rk)) .* tauk   ;
end



function iNDmix = LocalMixtureInternalEnergy(del,tau)
    
    % Get the saturation state
    [~,delL,delG] = SaturationStateGivenTauRRND(tau);
    
    % Calculate the mixture internal energy
    x     = QualityFromDensity(del,delL,delG);
    iLGND = InternalEnergyOneRND([delL;delG],[tau;tau]);
    iLND  = iLGND(1:end/2)      ;
    iGND  = iLGND(end/2+1:end)  ;
    iNDmix = iLND + x.*(iGND - iLND);
    
end





