function varargout = SaturationStateGivenMixedRhoIRRND(delMix,iMixND,varargin)
    
    %   Input patterns:
    %       o SaturationStateGivenRhoImixRRND(delMix,iMixND,tauGuess)
    %       o SaturationStateGivenRhoImixRRND(delMix,iMixND,tauGuess,tauSat)
    %       o SaturationStateGivenRhoImixRRND(delMix,iMixND,tauGuess,tauSat,iSatND)
    
	iSat = InternalEnergyOneR(delMix,tauSat);
    
    % Get initial guess for the iterative solution
    if (nargin < 3) || isempty(Tguess)
        %   The guess function performs a five-point interpolation on the closed
        %   temperature interval [Tt,Tsat].
        Nnodes = 2;
        tauGuess = GetTauGuess(delMix,iMix,delL,delG,tauSat,Nnodes);
    else
        tauGuess = Tc ./ Tguess;  % It may be user-supplied.
    end

    % Iteration Setup
    Tolerance = 1E-12;
    IterMax   = 1E3;
    UpdateFun = @(tau,Mask) Updater(tau,Mask,delMix,iMix,tauSat,iSat);
    Guess     = [tauGuess,0*tauGuess,0*tauGuess];
    
    % Newton iterations
    Solution = NewtonUpdater(UpdateFun,Guess,Tolerance,IterMax);
    
    % Pull the solution values
    tau  = Solution(:,1);
    delL = Solution(:,2);
    delG = Solution(:,3);
    
    % Compute output values
    Tsat = Tc   ./tau;
    rhol = delL .*rhoc;
    rhog = delG .*rhoc;
    Psat = PressureOneR(delL,tau);
    x    = QualityFromDensity(delMix,delL,delG);
    
    varargout = {Psat,delL,delG,Tsat,x};
    
end
