function [Psat,tauSat,delL,delG] = SaturationStateGivenDelta(delKn,UniqueMask)
    
    if (nargin < 2) || isempty(UniqueMask)
        [delKn,~,UniqueMask] = unique(delKn);
    end
    
    delL   = delKn*0;
    delG   = delKn*0;
    tauSat = delKn*0;
    
    LiquidKnown = delKn > 1;
    GasKnown    = delKn < 1;
    
    delLKn = delKn(LiquidKnown);
    delGKn = delKn(GasKnown);
    
    I      = 1:length(delKn);
    IpushG = I(LiquidKnown);
    IpushL = I(GasKnown   );
    
    Tolerance   = 1E-10;
    IterMax     = 1E3;
    
    % Solve the system using the liquid density as a known
    if any(LiquidKnown)
        [~,delGUn,tauGUn] = EstimateDelLDelGTauFromDel(delLKn);
        
        Calculate = delLKn >= DoNoIterationValueReducedDensityLiquid();
        
        [delGSol,tauGSol]   = SolveGasSystem(delLKn(Calculate),...
                                             delGUn(Calculate),...
                                             tauGUn(Calculate),...
                                             Tolerance,IterMax);
        delGUn(Calculate) = delGSol;
        tauGUn(Calculate) = tauGSol;
        
        delL  (IpushG) = delLKn;
        delG  (IpushG) = delGUn;
        tauSat(IpushG) = tauGUn;
    end
    
    % Solve the system using the gas density as a known
    if any(GasKnown)
        [delLUn,~,tauLUn] = EstimateDelLDelGTauFromDel(delGKn);
        
        Calculate = delGKn <= DoNoIterationValueReducedDensityGas() ;
        [delLSol,tauLSol] = SolveLiquidSystem(delLUn(Calculate),...
                                              delGKn(Calculate),...
                                              tauLUn(Calculate),...
                                              Tolerance,IterMax);

        delLUn(Calculate) = delLSol;
        tauLUn(Calculate) = tauLSol;

        delG  (IpushL) = delGKn;
        delL  (IpushL) = delLUn;
        tauSat(IpushL) = tauLUn;
    end
    
    % Assign the values of a critical system
    IsCritical         = not(GasKnown) & not(LiquidKnown);
    delL  (IsCritical) = 1;
    delG  (IsCritical) = 1;             
    tauSat(IsCritical) = 1;
    
%     % Produce output
%     Tsat =     Tc     ./ tauSat;
%     rhol = delL .* rhoc;
%     rhog = delG .* rhoc;
    Psat = PressureOneR(delL,tauSat);
    
    Psat   = Psat  (UniqueMask);
    tauSat = tauSat(UniqueMask);
    delL   = delL  (UniqueMask);
    delG   = delG  (UniqueMask);
end

function [delG,tau] = SolveGasSystem(delL,delG,tau,Tolerance,IterMax)
    
    % Get the first updates and residual values
    [RNorm,delGUpdate,tauUpdate] = UpdateGasSystem(delL,delG,tau);
    
    % First convergence check; simple L1-absolute check against the residuals.
    NotConvergedVec = RNorm > Tolerance   ;
    NotConverged    = any(NotConvergedVec);
    
    % Initialize iteration count
    Iter   = 0;
    
    % Build the vector of indices used for the push/contract scheme in the loop
    Icon   = 1:length(tau);
    Icon   = Icon(NotConvergedVec);
    
    % Contract already-converged vectors
    delLk = delL(NotConvergedVec);
    delGk = delG(NotConvergedVec);
    tauk  = tau (NotConvergedVec);
    
    while NotConverged && (Iter < IterMax)
        
        % Update the iteration variables
        delGk = delGk - delGUpdate(NotConvergedVec) ;
         tauk =  tauk -  tauUpdate(NotConvergedVec) ;
        
        % Calculate the next update
        [RNorm,delGUpdate,tauUpdate] = UpdateGasSystem(delLk,delGk,tauk);
        
        % Convergence check
        NotConvergedVec = RNorm > Tolerance     ;
        ConvergedVec    = not(NotConvergedVec)  ;
        NotConverged    = any(NotConvergedVec)  ;
                
        % Push the converged quantities to the output vectors
        Ipush       = Icon (ConvergedVec);
        delG(Ipush) = delGk(ConvergedVec);
        tau(Ipush)  =  tauk(ConvergedVec);
        
        % Contract the vectors used in calculations
        Icon   = Icon (NotConvergedVec);
        delLk  = delLk(NotConvergedVec);
        delGk  = delGk(NotConvergedVec);
        tauk   =  tauk(NotConvergedVec);
        
        Iter = Iter + 1;
    end
end

function [delL,tau] = SolveLiquidSystem(delL,delG,tau,Tolerance,IterMax)
    
    % Get the first updates and residual values
    [RNorm,delLUpdate,tauUpdate] = UpdateLiquidSystem(delL,delG,tau);
    
    % First convergence check; simple L1-absolute check against the residuals.
    NotConvergedVec = RNorm > Tolerance     ;
    NotConverged    = any(NotConvergedVec)  ;
    
    % Initialize iteration count
    Iter   = 0;
    
    % Build the vector of indices used for the push/contract scheme in the loop
    Icon   = 1:length(tau);
    Icon   = Icon(NotConvergedVec);
    
    % Contract already-converged vectors
    delLk = delL(NotConvergedVec);
    delGk = delG(NotConvergedVec);
    tauk  = tau (NotConvergedVec);
    
    while NotConverged && (Iter < IterMax)
        
        % Update the iteration variables
        delLk = delLk - delLUpdate(NotConvergedVec) ;
         tauk =  tauk -  tauUpdate(NotConvergedVec) ;
        
        % Calculate the next update
        [RNorm,delLUpdate,tauUpdate] = UpdateLiquidSystem(delLk,delGk,tauk);
        
        % Convergence check
        NotConvergedVec = RNorm > Tolerance     ;
        ConvergedVec    = not(NotConvergedVec)  ;
        NotConverged    = any(NotConvergedVec)  ;
        
        % Push the converged quantities to the output vectors
        Ipush       = Icon (ConvergedVec);
        delL(Ipush) = delLk(ConvergedVec);
        tau(Ipush)  =  tauk(ConvergedVec);
        
        % Contract the vectors used in calculations
        Icon   = Icon (NotConvergedVec);
        delLk  = delLk(NotConvergedVec);
        delGk  = delGk(NotConvergedVec);
        tauk   =  tauk(NotConvergedVec);
        
        Iter = Iter + 1;
    end
end

function [RNorm,ddelG,dtau] = UpdateGasSystem(delL,delG,tau)
    PhiR    = @(del,tau) HelmholtzResidual   (del,tau);
    PhiR_d  = @(del,tau) HelmholtzResidual_d (del,tau);
    PhiR_t  = @(del,tau) HelmholtzResidual_t (del,tau);
    PhiR_dd = @(del,tau) HelmholtzResidual_dd(del,tau);
    PhiR_dt = @(del,tau) HelmholtzResidual_dt(del,tau);
    
%     Pnodim = @(del,tau) 1 + del.*PhiR_d(del,tau);
    
    % Psig is the (essentiall) the Gibb's Energy Equality the is eliminated from
    % the system to shrink the solution space.  Here are the functions and
    % derivatives that are needed in the residuals below.
    %
    % Psig = Psig1 * Psig2 to more easily maintain/check the handles.
    Psig1   = PhiR(delL,tau) - PhiR(delG,tau) + log(delL./delG);
    Psig1_d = -1./delG - PhiR_d(delG,tau)                      ;
    Psig1_t = PhiR_t(delL,tau) - PhiR_t(delG,tau)              ;
    
    Psig2   = (delL.*delG)./(delL - delG)                      ;
    Psig2_d = ( delL./(delG - delL) ).^2                       ;
    % Psig2_t = 0;
    
    Psig    = Psig1 .* Psig2                        ;
    Psig_d  = Psig1_d .* Psig2 + Psig2_d .* Psig1   ;
    Psig_t  = Psig1_t .* Psig2                      ;
    
    
    %  Residuals from explicitly solving the 2x2 Newton update system.
    %
    R1   = Psig - delL.*(1 + delL .* PhiR_d(delL,tau));
    R2   = Psig - delG.*(1 + delG .* PhiR_d(delG,tau));
    R1_d = Psig_d;
    R1_t = Psig_t - delL.^2 .* PhiR_dt(delL,tau)   ;
    R2_d = Psig_d - (1 + 2*delG.*PhiR_d(delG,tau) + delG.^2.*PhiR_dd(delG,tau));
    R2_t = Psig_t - delG.^2 .* PhiR_dt(delG,tau)   ;
    DetR =  R2_t .* R1_d - R1_t .* R2_d ;
    
    
    % Final update formula
    ddelG   = (R1 .* R2_t - R2 .* R1_t) ./ DetR  ;
    dtau    = (R2 .* R1_d - R1 .* R2_d) ./ DetR  ;
    RNorm   = abs(R1) + abs(R2);
end

function [RNorm,ddelL,dtau] = UpdateLiquidSystem(delL,delG,tau)
    PhiR    = @(del,tau) HelmholtzResidual   (del,tau);
    PhiR_d  = @(del,tau) HelmholtzResidual_d (del,tau);
    PhiR_t  = @(del,tau) HelmholtzResidual_t (del,tau);
    PhiR_dd = @(del,tau) HelmholtzResidual_dd(del,tau);
    PhiR_dt = @(del,tau) HelmholtzResidual_dt(del,tau);
    
   
    % Psig is the (essentiall) the Gibb's Energy Equality the is eliminated from
    % the system to shrink the solution space.  Here are the functions and
    % derivatives that are needed in the residuals below.
    %
    % Psig = Psig1 * Psig2 to more easily maintain/check the handles.
    Psig1   = PhiR(delL,tau) - PhiR(delG,tau) + log(delL./delG) ;
    Psig1_d = 1./delL + PhiR_d(delL,tau)                        ;
    Psig1_t = PhiR_t(delL,tau) - PhiR_t(delG,tau)               ;
    
    Psig2   = (delL.*delG)./(delL - delG)   ;
    Psig2_d = -( delG./(delG - delL) ).^2   ;
    % Psig2_t = 0;
    
    Psig    = Psig1   .* Psig2                      ;
    Psig_d  = Psig1_d .* Psig2 + Psig2_d .* Psig1   ;
    Psig_t  = Psig1_t .* Psig2                      ;
    
    
    %  Residuals from explicitly solving the 2x2 Newton update system.
    %
    R1   = Psig  - delL .* (1 + delL.*PhiR_d(delL,tau));
    R2   = Psig  - delG .* (1 + delG.*PhiR_d(delG,tau));
    R1_d = Psig_d - (1 + 2*delL.*PhiR_d(delL,tau) + delL.^2.*PhiR_dd(delL,tau));
    R1_t = Psig_t - delL.^2 .* PhiR_dt(delL,tau)    ;
    R2_d = Psig_d                                   ;
    R2_t = Psig_t - delG.^2 .* PhiR_dt(delG,tau)    ;
    DetR = R1_d .* R2_t - R1_t.*R2_d                ;
    
    
    % Final update formula
    ddelL =  (R1 .* R2_t - R2 .* R1_t) ./ DetR ;
    dtau  =  (R2 .* R1_d - R1 .* R2_d) ./ DetR ;
    RNorm = abs(R1) + abs(R2);
end

% function tau = TauGuessGasSide(del)
%     c = [+3.104332488E-004,...
%          +3.975447616E-003,...
%          +6.202202323E-003,...
%          +3.788570083E-002,...
%          -3.575927611E-001,...
%          +1.441183511E+000];
%         
%     mu = [+7.029298622E-001,+1.458251049E-001];
%     
%     tau =  HornersMethod(del.^(+8.341717E-2),c,mu);    
% end
% 
% function tau = TauGuessLiquidSide(del)
%     c = [+3.226566684E-002,...
%          -1.555966390E-002,...
%          -2.936292390E-002,...
%          +2.788300063E-002,...
%          +3.043589498E-001,...
%          +1.455237966E+000];
%         
%     mu = [+2.687134584E+003,+2.074491879E+003];
%     
%     tau =  HornersMethod(del.^7.714444,c,mu);    
% end
% 
% function dell = LiquidGuess(del)
%     c = [-4.623661465E-003,...
%          -2.611856179E-002,...
%          -6.966175395E-002,...
%          -1.702069908E-001,...
%          -3.444782271E-001,...
%          +2.803200913E+000];
%         
%     mu = [+9.649749715E-001,+2.016717391E-002];
%     
%     dell =  HornersMethod(del.^7.955455E-003,c,mu);    
% end
% 
% function delg = GasGuess(del)
%     c = [-2.274370237E-004,...
%          -8.065695223E-004,...
%          -5.876864143E-003,...
%          +4.038671626E-002,...
%          -6.462268409E-002,...
%          +3.180756392E-002];
%         
%     mu = [+3.703732222E+000,+7.739607151E-001];
%     
%     delg =  HornersMethod(del.^1.351414,c,mu);   
% end














