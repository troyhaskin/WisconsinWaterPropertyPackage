function varargout = SaturationStateGivenMixedRhoIRRND(delMix,iNDmix,varargin)
    
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
    
    
    Sizer = zeros(size(delMix));
    
    
    % Calculate the mixture internal energy on the triple line
    [delLt,delGt] = TriplePointDensitiesR()                 ;
    [iLNDt,iGNDt] = TriplePointInternalEnergiesND()         ;
    taut          = TriplePointTemperatureR() + Sizer       ;
    xt            = QualityFromDensity(delMix,delLt,delGt)  ;
    iNDt          = iLNDt + xt.*(iGNDt - iLNDt)             ;
    
    
    % Default/optional arguments
    OptionalArguments = {true,[],[],[]}             ;   % Allocate with defaults
    nOptionalIn       = length(varargin)            ;   % Number of optionals passed
    OptionalArguments(1:nOptionalIn) = varargin     ;   % Overwrite defaults
    [PhaseCheck,tauSat,iSatND,tauGuess] = OptionalArguments{:};
    
    % Default phase check
    if isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    % If no tauSat is given for the delMix, calculate it.
    if isempty(tauSat)
        [~,tauSat,~,~] = SaturationStateGivenDelta(delMix);
    end
    
    % If no iSatND is given for the delMix, calculate it.
    if isempty(iSatND)
        iSatND = InternalEnergyOneRND(delMix,tauSat);
    end
    
    % If no tauGuess is given, calculate it from a secant line.
    if isempty(tauGuess)
        tauGuess = taut + (iNDmix - iNDt)./(iSatND - iNDt) .* (tauSat - taut);
    end
    
    
    % Iteration Setup
    Tolerance     = 1E-12 ;
    IterMax       = 1E2   ;
    NotDone       = true  ;
    Iter          = 0     ;
    tauSol        = Sizer ;
    iNDwork       = iNDmix;
    delWork       = delMix;
    InotConverged = 1:length(Sizer);
    
    % Define the first three iterates and the nex
    taukm2 = tauSat     ;
    taukm1 = taut       ;
    tauk   = tauGuess   ;
    
    % Define the first three residuals
    Rkm2 = iSatND - iNDmix ;
    Rkm1 = iNDt   - iNDmix ;
    Rk   = LocalMixtureInternalEnergy(delMix,tauGuess) - iNDmix;

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
    
    [Psat,delL,delG] = SaturationStateGivenTausat(tauSol);
    x = QualityFromDensity(delMix,delL,delG);
    
    varargout = {Psat,tauSol,delL,delG,x};
end

function iNDmix = LocalMixtureInternalEnergy(del,tau)
    
    % Get the saturation state
    [~,delL,delG] = SaturationStateGivenTausat(tau);
    
    % Calculate the mixture internal energy
    x    = QualityFromDensity(del,delL,delG);
    iLND = InternalEnergyOneRND(delL,tau);
    iGND = InternalEnergyOneRND(delG,tau);
    iNDmix = iLND + x.*(iGND - iLND);
    
end





