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
    Tolerance   = 1E-12 ;
    IterMax     = 1E2   ;
    NotDone     = true  ;
    Iter        = 0     ;
    tauSol      = Sizer ;
%     taukm1      = Sizer ;
%     iNDkm1      = Sizer ;
    iNDwork     = iNDmix;
    delWork     = delMix;
% 
%     CloserToTauSat = (tauGuess - tauSat) < (taut - tauGuess);
%     CloserToTauT   = not(CloserToTauSat);
%     
%     taukm1(CloserToTauSat) = tauSat(CloserToTauSat) ;
%     iNDkm1(CloserToTauSat) = iSatND(CloserToTauSat) ;
%     taukm1(CloserToTauT)   = taut  (CloserToTauT  ) ;
%     iNDkm1(CloserToTauT)   = iNDt  (CloserToTauT  ) ;

    taukm2 = tauSat ;
    taukm1 = taut   ;

    iNDkm2 = iSatND - iNDmix ;
    iNDkm1 = iNDt   - iNDmix ;
    
    tauk   = tauGuess   ;
    iNDk   = LocalMixtureInternalEnergy(delMix,tauGuess) - iNDmix;

    while NotDone
%         taukp1 = tauk - (tauk - taukm1) .* (iNDk - iNDwork)./(iNDk - iNDkm1);

        taukp1 = iNDkm1.*iNDk  ./((iNDkm2-iNDkm1).*(iNDkm2-iNDk)) .* taukm2 + ...
                 iNDkm2.*iNDk  ./((iNDkm1-iNDkm2).*(iNDkm1-iNDk)) .* taukm1 + ...
                 iNDkm2.*iNDkm1./((iNDkm2-iNDk  ).*(iNDkm1-iNDk)) .* tauk   ;

        iNDkp1 = LocalMixtureInternalEnergy(delWork,taukp1) - iNDwork;

        Converged    = (abs(taukp1-tauk) < Tolerance) | (abs(iNDkp1-iNDk) < Tolerance);
        NotConverged = not(Converged);

        % Push converged values
        tauSol(Converged) = taukp1(Converged);

        % Contract all other vectors
        taukm2  = taukm1 (NotConverged);
        taukm1  = tauk   (NotConverged);
        tauk    = taukp1 (NotConverged);
        iNDkm2  = iNDkm1 (NotConverged);
        iNDkm1  = iNDk   (NotConverged);
        iNDk    = iNDkp1 (NotConverged);
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





