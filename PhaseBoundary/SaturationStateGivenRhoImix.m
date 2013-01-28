function [Psat,Tsat,rhol,rhog,x] = SaturationStateGivenRhoImix(rhoMix,iMix,Tguess,Tsat)
    
    % Reference values
    rhoc  = CriticalDensity();
    Tc    = CriticalTemperature();
    
    % Deal with vector sizes (mostly scalar vs. vector)
    [rhoMix,iMix]   = BalanceSizes(rhoMix,iMix) ;
    
    % Reduced quantities
    delMix = rhoMix  / rhoc ;

    % Get the saturation state from the reduced density
    if (nargin < 4) || isempty(Tsat)
        [~,tauSat,delL,delG] = SaturationStateGivenDelta(delMix) ; % Get tauSat from delMix
    else
        tauSat = Tc ./Tsat;  % It may be user-supplied.
        [~,delL,delG] = SaturationStateGivenTausat(tauSat) ;
    end
    
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
    
end

function [dSys,RNorm] = Updater(Sys,Mask,deltaMix,iMix,tauSat,iSat)
    
    tau  = Sys(:,1);
    
    istar = DimensioningInternalEnergy();
    
    % Pull only non-converged values
    delGiven       = deltaMix(Mask) ;
    IntEnergyGiven = iMix(Mask)     ;
    tauSatGiven    = tauSat(Mask)   ;
    iSatGiven      = iSat(Mask)     ;
    
    % Get the saturation state
    [Psat,delL,delG] = SaturationStateGivenTausat(tau);
    
    % Calculate the mixture internal energy
    x  = QualityFromDensity(delGiven,delL,delG);
    il = InternalEnergyOneR(delL,tau);
    ig = InternalEnergyOneR(delG,tau);
    IntEnergy = il + x.*(ig - il);
    
    % Total derivative of mixture internal energy w.r.t. tau
 	IntEnergy_tau = InternalEnergyR_tauSat(tau,delGiven,Psat,delL,delG,il,ig); % This is technically wrong, but it gets the job done faster.
%     IntEnergy_tau = InternalEnergyR_tauN(delGiven,tau); 
    
    % Form the residuals
    Residual    = (IntEnergy - IntEnergyGiven)/istar  ;
    dResidual   = IntEnergy_tau / istar               ;
    
    % Calculate the update for tau and it's norm
    dtau  = Residual ./ dResidual   ;
    RNorm = abs(Residual)           ;
    
    
    % Guard against flying outside the vapor dome (since we know we're in there)
    tNewton   = tau - dtau;
    OutsideVP = tNewton < tauSatGiven;
    if any(OutsideVP)
        m               = (IntEnergy(OutsideVP) - IntEnergyGiven(OutsideVP))    ./...
                          (IntEnergy(OutsideVP) - iSatGiven(OutsideVP))         ;
        dtau(OutsideVP) = m.*(tau(OutsideVP) - tauSatGiven(OutsideVP))          ;
    end
    
    
    % The ONLY part of the system being updated is tau, but in order to pull
    % the reduced densities so as not to recalculate them after the Newton 
    % update is finished, this fake update is used that maintains the densities
    % calculated from the current tau iterant.
    if (Sys(:,2) ~= 0)
        ddel = Sys(:,2:3) - [delL,delG];
        dSys = [dtau,ddel];
    else
        dSys = [dtau,-delL,-delG];
    end
    
end


function tauGuess = GetTauGuess(delMix,iMix,delL,delG,tauSat,Nnodes)
    
    %   Triple point information
	%       By convention for water, ilt === 0, so it is simply ignored.
    taut          = CriticalTemperature() / TriplePointTemperature();
    [delLt,delGt] = TriplePointReducedDensities();  % Saturated densities for quality
    [~,igt]       = TriplePointInternalEnergies();  % Saturated internal energies
    
    % Triple line internal energy
    xt = QualityFromDensity(delMix,delLt,delGt);    % Quality
    it = xt .* igt;                                 % Mixture Triple line internal energy
    
    % Saturation internal energy
    isat = InternalEnergyMixtureR(delMix,delL,delG,tauSat);
    
    
    % Normalized interpolation nodes
    alpha = ChebyshevSpace(0,1,Nnodes); % Weights for interior points
%     alpha = LinearSpace(0,1,Nnodes); % Weights for interior points
    
    % Interpolation arrays 
    % tau nodes 
    tauInterp = zeros(length(tauSat),Nnodes);
    tauInterp(:,1)   = tauSat;
    tauInterp(:,end) = CriticalTemperature()/TriplePointTemperature(); % Triple point tau
    
    % internal energy nodes
    iInterp = zeros(length(tauSat),Nnodes);
    iInterp(:,1)   = isat;
    iInterp(:,end) = it;
    
    % Fill in the interior data
    for k = 2:Nnodes-1
        % Interior tau interpolation node
        tauMid = tauSat + alpha(k).*(taut - tauSat);
        
        % Saturation state
        [~,delL,delG]  = SaturationStateGivenTausat(tauMid);
       
        % Mixture internal energy at this node
        iInterp(:,k)   = InternalEnergyMixtureR(delMix,delL,delG,tauMid);
        
        % Push the tau values
        tauInterp(:,k) = tauMid;
    end

    % Lagrange interpolation of the node data above
    istar = DimensioningInternalEnergy();
    tauGuess = LagrangeInterpolation(iMix/istar,iInterp/istar,tauInterp);
    
end


