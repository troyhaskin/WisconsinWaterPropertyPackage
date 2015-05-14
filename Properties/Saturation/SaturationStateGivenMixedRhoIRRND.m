function [Pnd,tauSol,delL,delG,x] = SaturationStateGivenMixedRhoIRRND(delta,iND,varargin)
    
    %   Input patterns:
    %       o SaturationStateGivenMixedRhoIRRND(delMix,PhaseCheck,iMixND)
    %       o SaturationStateGivenMixedRhoIRRND(delMix,PhaseCheck,iMixND,tauSat,iSatND)
    %       o SaturationStateGivenMixedRhoIRRND(delMix,PhaseCheck,iMixND,tauSat,iSatND,tauGuess)
    %
    %   For this function, properties explicitly coded as <S/sat> refer to purely saturated
    %   (x = 0 or 1 only) properties on the delMix isochor (delSat = delMix).  While all
    %   tau greater than tauSat are also in saturation, tauSat is the smallest tau value
    %   at which a two-phase state can exist (while the largest is the triple-line state);
    %   therefore, the <S/sat> values allow the secant method to bounded/guarded.
    %
    
    
    Sizer = zeros(size(delta));
    
    
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
    
    % ============================================ %
    %               Exact Phase Check              %
    % ============================================ %
    
    %   non-DMVS region
    lessThanTriple = delta < delLt;
    [~,tauHi(lessThanTriple),~,~] = SaturationStateGivenDeltaRRND(delta(lessThanTriple));
    tauLo(lessThanTriple)         = taut                                                ;
    
    %   DMVS region
    inDMVS                   = not(lessThanTriple);
    deltaDMVS                = delta(inDMVS);
    tauBotTop                = estimateTauDMVSRegion(deltaDMVS);
    [~,tauNearTriple,~,~]    = SaturationStateGivenDeltaRRND([deltaDMVS;deltaDMVS],tauBotTop);
    tauLo(inDMVS) = tauNearTriple(   1   :end/2);
    tauHi(inDMVS) = tauNearTriple(end/2+1: end );

    %   Check upper bound
    iHi          = InternalEnergyOneRND(delta,tauLo) ;
    notSaturated = iND > iHi                         ;
    
    %   Check lower bound
    iLo = InternalEnergyOneRND(deltaDMVS,tauLo(inDMVS)) ;
    notSaturated(inDMVS) = notSaturated(inDMVS) | (iND(inDMVS) < iLo(inDMVS));
    
    
    % If no iSatND is given for the delMix, calculate it.
    if isempty(iSatND)
        iSatND = InternalEnergyOneRND(delta,tauSat);
    end
    
    % If no tauGuess is given, calculate it from a secant line.
    if isempty(tauGuess)
        tauGuess = taut + (iND - iNDt)./(iSatND - iNDt) .* (tauSat - taut);
    end
    
    
    % Iteration Setup
    Tolerance     = 1E-12 ;
    IterMax       = 1E2   ;
    NotDone       = true  ;
    Iter          = 0     ;
    tauSol        = Sizer ;
    iNDwork       = iND;
    delWork       = delta;
    InotConverged = 1:length(Sizer);
    
    % Define the first three iterates and the nex
    taukm2 = tauSat     ;
    taukm1 = taut       ;
    tauk   = tauGuess   ;
    
    % Define the first three residuals
    Rkm2 = iSatND - iND ;
    Rkm1 = iNDt   - iND ;
    Rk   = LocalMixtureInternalEnergy(delta,tauGuess) - iND;

% ======================================================================= %
%                            Solution Block                               %
% ======================================================================= %
    % Check for zeroth-guess convergence
    Converged    = abs(Rk) < Tolerance;
    NotConverged = not(Converged);
    
    % Update converged index array
    Ipush         = InotConverged(Converged)    ;
    InotConverged = InotConverged(NotConverged) ;
    tauSol(Ipush) = tauk(Converged)             ;

    % Contract all other vectors
    taukm2  = taukm2 (NotConverged);
    taukm1  = taukm1 (NotConverged);
    tauk    = tauk   (NotConverged);
    Rkm2    = Rkm2   (NotConverged);
    Rkm1    = Rkm1   (NotConverged);
    Rk      = Rk     (NotConverged);
    iNDwork = iNDwork(NotConverged);
    delWork = delWork(NotConverged);
    
    % Enter Iterative loop
    while NotDone
        
        taukp1 = Rkm1.*Rk  ./((Rkm2-Rkm1).*(Rkm2-Rk)) .* taukm2 + ...
                 Rkm2.*Rk  ./((Rkm1-Rkm2).*(Rkm1-Rk)) .* taukm1 + ...
                 Rkm2.*Rkm1./((Rkm2-Rk  ).*(Rkm1-Rk)) .* tauk   ;
        
        Rkp1 = LocalMixtureInternalEnergy(delWork,taukp1) - iNDwork;
        
        Converged    = (abs(taukp1-tauk) < Tolerance) | (abs(Rkp1) < Tolerance);
        NotConverged = not(Converged);
        
        % Update index array
        Ipush         = InotConverged(Converged);
        InotConverged = InotConverged(NotConverged);
        
        % Push converged values
        tauSol(Ipush) = taukp1(Converged);
        
        % Contract all other vectors
        taukm2  = taukm1 (NotConverged);
        taukm1  = tauk   (NotConverged);
        tauk    = taukp1 (NotConverged);
        Rkm2    = Rkm1   (NotConverged);
        Rkm1    = Rk     (NotConverged);
        Rk      = Rkp1   (NotConverged);
        iNDwork = iNDwork(NotConverged);
        delWork = delWork(NotConverged);
        
        % Loop break check
        Iter = Iter + 1;
        NotDone = any(NotConverged) && (Iter < IterMax);
    end
    
    [Pnd,delL,delG] = SaturationStateGivenTauRRND(tauSol);
    x = QualityFromDensity(delta,delL,delG);

end

function iNDmix = LocalMixtureInternalEnergy(del,tau)
    
    % Get the saturation state
    [~,delL,delG] = SaturationStateGivenTauRRND(tau);
    
    % Calculate the mixture internal energy
    x    = QualityFromDensity(del,delL,delG);
    iLND = InternalEnergyOneRND(delL,tau);
    iGND = InternalEnergyOneRND(delG,tau);
    iNDmix = iLND + x.*(iGND - iLND);
    
end





