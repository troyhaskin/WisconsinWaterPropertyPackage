function [Psat,Tsat,rhoL,rhoG,x] = SaturationStateGivenMixedRhoI(rhoMix,iMix,varargin)
    
    % Dimensionalizers
    [~,rhoc,Tc] = Nondimensionalizers();
    istar       = DimensioningInternalEnergy();
    
    delMix = rhoMix / rhoc  ;
    iNDmix = iMix   / istar ;

    % Default/optional arguments
	OptionalArguments = {true,[],[],[]}             ;   % Allocate with defaults
    nOptionalIn       = length(varargin)            ;   % Number of optionals passed
    OptionalArguments(1:nOptionalIn) = varargin     ;   % Overwrite defaults
    PhaseCheck = OptionalArguments{1};
    tauSat     = Tc ./ OptionalArguments{2};
    tauGuess   = Tc ./ OptionalArguments{3};
    iSatND     = OptionalArguments{3} / istar;


    % Call lower level function
    [Psat,tauSol,delL,delG,x] = SaturationStateGivenMixedRhoIRRND(...
                                        delMix,iNDmix,PhaseCheck,tauSat,iSatND,tauGuess);

    % Redimensionalize
    Tsat = Tc   ./ tauSol   ;
    rhoL = delL  * rhoc     ;
    rhoG = delG  * rhoc     ;

end

