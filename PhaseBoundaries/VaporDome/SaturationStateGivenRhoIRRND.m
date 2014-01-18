function varargout = SaturationStateGivenRhoIRRND(del,iND,varargin)
    
    %   Input patterns:
    %       o SaturationStateGivenRhoImixRRND(delMix,iMixND,tauGuess)
    %       o SaturationStateGivenRhoImixRRND(delMix,iMixND,tauGuess,tauSat)
    %       o SaturationStateGivenRhoImixRRND(delMix,iMixND,tauGuess,tauSat,iSatND)
    
    [~,tauGuess,~,~] = SaturationStateGivenDelta(del);
    
    % Iteration Setup
    UpdateFun = @(tau,Mask) Updater(tau,Mask,del,iND)   ;
    Guess     = tauGuess    ;
    IterMax   = 1E3         ;
    Tolerance = 1E-12       ;
    
    % Newton iterations
    tau = NewtonUpdater(UpdateFun,Guess,Tolerance,IterMax);
    
    % Compute output values
    [Pnd,delL,delG] = SaturationStateGivenTauRRND(tau);
    x    = QualityFromDensity(del,delL,delG);
    
    % Output
    varargout = {Pnd,delL,delG,tau,x};
    
end


function [dtau,Rnorm] = Updater(tau,Mask,del,iND)
    
    
    iND = iND(Mask);
    del = del(Mask);
    
    [Pnd,delL,delG] = SaturationStateGivenTauRRND(tau);
    
    iL  = InternalEnergyOneRND(delL,tau);
    iG  = InternalEnergyOneRND(delG,tau);
    x   = delG./del .* (del - delL)./(delG - delL);
    
    % Phase determination
    OnePhase = (x >= 0.9999) | (x <= 0.0001)  ;
    TwoPhase = not(OnePhase)        ;
    
    % Allocation
    R  = x*0    ;
    dR = R      ;
    
    if any(OnePhase)
        R1  = InternalEnergyOneRND(del(OnePhase),tau(OnePhase)) - iND(OnePhase);
        dR1 = InternalEnergyOneRND_tau(del(OnePhase),tau(OnePhase));
        
        % Assignment
        R (OnePhase) = R1   ;
        dR(OnePhase) = dR1  ;
    end
    
    if any(TwoPhase)
        [R2,dR2] = TwoPhaseUpdater(Pnd (TwoPhase),tau (TwoPhase)  ,...
            delL(TwoPhase),delG(TwoPhase)  ,....
            iL  (TwoPhase),iG  (TwoPhase)  ,...
            x(TwoPhase),del(TwoPhase),iND(TwoPhase));

        % Assignment
        R (TwoPhase) = R2   ;
        dR(TwoPhase) = dR2  ;
    end

    % Newton step
    dtau = R./dR;
    
    % L_1 residual
    Rnorm = abs(R);
    
end


function [R,dR] = TwoPhaseUpdater(Pnd,tau,delL,delG,iL,iG,x,del,iND)
    
    % Latent Internal Energy
    iLG = iG - iL;
    
    % Residual
    R   = iL + x.* iLG - iND;
    
    % Derivative of iL w.r.t. tau
    Pnd_tau  = ClausiusClapeyronRRND(Pnd,tau,delL,delG,iL,iG);
    PhiR_d   = HelmholtzResidual_d (delL,tau);
    PhiR_dd  = HelmholtzResidual_dd(delL,tau);
    PhiR_dt  = HelmholtzResidual_dt(delL,tau);
    delL_tau = (delL + tau.^2 .* Pnd_tau + delL.^2 .* (PhiR_d - tau .* PhiR_dt)) ./ ...
        (tau .* (1 + 2 * delL .* PhiR_d + delL.^2 .* PhiR_dd ));
    diL = InternalEnergyOneRND_delta(delL,tau) .* delL_tau + InternalEnergyOneRND_tau(delL,tau);
    
    % Derivative of iG w.r.t. tau
    PhiR_d   = HelmholtzResidual_d (delG,tau);
    PhiR_dd  = HelmholtzResidual_dd(delG,tau);
    PhiR_dt  = HelmholtzResidual_dt(delG,tau);
    delG_tau = (delG + tau.^2 .* Pnd_tau + delG.^2 .* (PhiR_d - tau .* PhiR_dt)) ./ ...
        (tau .* (1 + 2 * delG .* PhiR_d + delG.^2 .* PhiR_dd ));
    diG = InternalEnergyOneRND_delta(delG,tau) .* delG_tau + InternalEnergyOneRND_tau(delG,tau);
    
    % Derivative of iLG w.r.t. tau
    diLG = diG - diL;
    
    % Derivative of x w.r.t. tau
    dx = ((delL-del).*delL.*delG_tau + (del-delG).*delG.*delL_tau)./(del.*(delG-delL).^2);
    
    % Derivative of residual w.r.t. tau
    dR = diL + x .* diLG + iLG .* dx;
    
end