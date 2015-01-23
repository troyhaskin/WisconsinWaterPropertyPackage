function [Pnd,delL,delG] = SaturationStateGivenTauRRND(tau,delL0,delG0,UniqueMask)
    
    % Uniqueness filter; if a UniqueMask is passed; tauSat is assumed to have unique entries.
    if (nargin < 4) || isempty(UniqueMask)
        [tau,~,UniqueMask] = unique(tau);
    end
    
    % Reduced liquid density guess value
    if (nargin < 2) || isempty(delL0)
        delL  = EstimateDelLFromTau(tau);
    else
        delL = delL0;
    end
    
    % Reduced gas density guess value
    if (nargin < 3) || isempty(delG0)
        delG  =  EstimateDelGFromTau(tau);
    else
        delG = delG0;
    end
    
    % Check tau for how close to critical it is
    NearTc     = tau < DoNoIterationTau()   ;
    AboveTc    = tau <   1                  ;   % Check if Tsat is above the critical point
    AboveTt    = tau <= TriplePointTau()    ;
    NotAboveTc = not(AboveTc)               ;
    BelowTt    = not(AboveTt)               ;
    
    
    % Iteration setup
    Tolerance = 1E-14                                   ; % Abolsute iteration tolerance
    IterMax   = 300                                     ; % Maximum iteration count
    Calculate = not(NearTc) & NotAboveTc & AboveTt      ; % Mask for temps. not close to the critical point
    Guess     = [delL(Calculate),delG(Calculate)]       ; % Starting values for the iteration
%     UpdateFun = @(x,Mask) Updater(x,Mask,tau(Calculate)); % Function used for updating the solution
    UpdateFun = @(x,Mask) UpdaterConjugateGradient(x,Mask,tau(Calculate)); % Function used for updating the solution
    
    % Solve the system
    if any(Calculate)
        xSol = NewtonUpdater(UpdateFun,Guess,Tolerance,IterMax);
        
        % Update the iterated values
        delL(Calculate) = xSol(:,1);
        delG(Calculate) = xSol(:,2);
        
    end
    
    % These values above for T above the vapor dome ensure that the fluid will
    % be flagged as single phase.
    delL(AboveTc | BelowTt) = 0.0;
    delG(AboveTc | BelowTt) = 0.0;
    
    % If exactly critical, assign exact values
    delL(tau == 1) = 1;
    delG(tau == 1) = 1;
    
    % Saturation pressure
    Pnd             = delG                                              ; % allocate
    Pnd(NotAboveTc) = PressureOneRND(delL(NotAboveTc),tau(NotAboveTc))  ;
    Pnd(tau == 1)   = CriticalPressureND()                              ;
    Pnd(AboveTc | BelowTt) = 0                                          ;
    
    % Expand vectors back to non-unique lengths
    Pnd  = Pnd (UniqueMask) ;
    delL = delL(UniqueMask) ;
    delG = delG(UniqueMask) ;
    
end


function [dx,RNorm] = Updater(x,Mask,tau0)
    
    N    = length(Mask);
    tau  = [tau0(Mask);tau0(Mask)];
    delL = x(:,1);
    delG = x(:,2);
    
    
    % Form Jacobian determinant and inverse
    [PhiR,PhiR_d,PhiR_dd] = HelmholtzResidualCombo__d_dd([delL;delG],tau);
    [PhiRL,PhiRG,PhiR_dL,PhiR_dG,PhiR_ddL,PhiR_ddG] = VectorChunk([PhiR;PhiR_d;PhiR_dd],N);
    
    Psig1    =  PhiRL - PhiRG + log(delL./delG) ;
    Psig1_dL =   PhiR_dL + 1./delL              ;
    Psig1_dG = -(PhiR_dG + 1./delG)             ;
    
    PdelL = 1 + delL.*PhiR_dL;
    PdelG = 1 + delG.*PhiR_dG;
    
    PdelL_d = PhiR_dL + delL.*PhiR_ddL;
    PdelG_d = PhiR_dG + delG.*PhiR_ddG;
    
    R1 = delG.*Psig1 - PdelL.*(delL-delG);
    R2 = delL.*Psig1 - PdelG.*(delL-delG);
    
    R1_dL = delG.*Psig1_dL + (delG-delL).*PdelL_d - PdelL;
    R1_dG = Psig1 + PdelL  + delG.*Psig1_dG;
    
    R2_dL = Psig1 - PdelG  + delL.*Psig1_dL;
    R2_dG = delL.*Psig1_dG + (delG-delL).*PdelG_d + PdelG;
    
    % Determinant
    DetJ = R1_dL .* R2_dG - R1_dG .*R2_dL;
    
    % Inverse Jacobian elements
    iJ11 =  R2_dG ./ DetJ;
    iJ12 = -R1_dG ./ DetJ;
    iJ21 = -R2_dL ./ DetJ;
    iJ22 =  R1_dL ./ DetJ;



    % Newton updates
    ddelL = iJ11 .* R1 + iJ12 .* R2;
    ddelG = iJ21 .* R1 + iJ22 .* R2;
    
    % Pack updates and calculate norm for the Newton updater
    dx    = [ddelL,ddelG];
    RNorm = abs(R1) + abs(R2);
    
    Show(max(RNorm));
end

% function [delL,delG] = UpdaterBroyden(delLk,delGk,tau)
%
%     % HFE handles
%     PhiR    = @(del) HelmholtzResidual   (del,tau);
%     PhiR_d  = @(del) HelmholtzResidual_d (del,tau);
%     PhiR_dd = @(del) HelmholtzResidual_dd(del,tau);
%
%
%     epsilon  = 1E-7;
%
%     % Dimensionless pressure
%     pi = @(delL,delG) (PhiR(delL) - PhiR(delG) + log(delL./delG))./(tau.*(1./delG-1./delL));
%
%     RL = @(delL,delG) tau .* pi(delL,delG) - delL .* (1 + delL .* PhiR_d(delL));
%     RG = @(delL,delG) tau .* pi(delL,delG) - delG .* (1 + delG .* PhiR_d(delG));
%
%     dRL = @(delL,delG) (RL(delL+epsilon,delG) - RL(delL-epsilon,delG))./(2*epsilon);
%     dRG = @(delL,delG) (RG(delL,delG+epsilon) - RG(delL,delG-epsilon))./(2*epsilon);
%
%     ddRL = @(delL,delG) (RL(delL+epsilon,delG) - 2*RL(delL,delG)+RL(delL-epsilon,delG))./epsilon^2;
%     ddRG = @(delL,delG) (RG(delL,delG+epsilon) - 2*RG(delL,delG)+RG(delL,delG-epsilon))./epsilon^2;
%
%     while true
%
%         delLkm1 = delLk;
%         delGkm1 = delGk;
%
%         ddelLk = 2 * RL(delLkm1,delGkm1) .* dRL(delLkm1,delGkm1) ./ ...
%             (RL(delLkm1,delGkm1) .* ddRL(delLkm1,delGkm1) - 2 * dRL(delLkm1,delGkm1).^2);
%
%         ddelGk = 2 * RG(delLkm1,delGkm1) .* dRG(delLkm1,delGkm1) ./ ...
%             (RG(delLkm1,delGkm1) .* ddRG(delLkm1,delGkm1) - 2 * dRG(delLkm1,delGkm1).^2);
%
%         R = abs(RL(delLkm1,delGkm1)) + abs(RG(delLkm1,delGkm1));
%         Show([delLkm1;delGkm1;abs(RL(delLkm1,delGkm1)) + abs(RG(delLkm1,delGkm1))]);
%
%         delLk = delLkm1 + 0.001 * ddelLk;
%         delGk = delGkm1 + 0.001 * ddelGk;
%
%         if R < 5E-8;
%             g = [];
%         end
%
%     end
%
% end

function [dx,RNorm] = UpdaterConjugateGradient(x,Mask,tau0)

    N    = length(Mask);
    tau  = [tau0(Mask);tau0(Mask)];
    delL = x(:,1);
    delG = x(:,2);
    
    
    % Form Jacobian determinant and inverse
    [PhiR,PhiR_d,PhiR_dd] = HelmholtzResidualCombo__d_dd([delL;delG],tau);
    [PhiRL,PhiRG,PhiR_dL,PhiR_dG,PhiR_ddL,PhiR_ddG] = VectorChunk([PhiR;PhiR_d;PhiR_dd],N);
    
    Psig1    =  PhiRL - PhiRG + log(delL./delG) ;
    Psig1_dL =   PhiR_dL + 1./delL              ;
    Psig1_dG = -(PhiR_dG + 1./delG)             ;
    
    PdelL = 1 + delL.*PhiR_dL;
    PdelG = 1 + delG.*PhiR_dG;
    
    PdelL_d = PhiR_dL + delL.*PhiR_ddL;
    PdelG_d = PhiR_dG + delG.*PhiR_ddG;
    
    R1 = delG.*Psig1 - PdelL.*(delL-delG);
    R2 = delL.*Psig1 - PdelG.*(delL-delG);
    
    R1_dL = delG.*Psig1_dL + (delG-delL).*PdelL_d - PdelL;
    R1_dG = Psig1 + PdelL  + delG.*Psig1_dG;
    
    R2_dL = Psig1 - PdelG  + delL.*Psig1_dL;
    R2_dG = delL.*Psig1_dG + (delG-delL).*PdelG_d + PdelG;
    

    % Conjugate gradient
    rk1 = R1;
    rk2 = R2;
    
    xk1 = 0 ;
    xk2 = 0 ;
    
    pk1 = R1;
    pk2 = R2;

    % First Step
    Ap1   = pk1 .* R1_dL + pk2 .* R1_dG;
    Ap2   = pk1 .* R2_dL + pk2 .* R2_dG;
    alpha = (rk1.^2 + rk2.^2)./(pk1 .* Ap1 + pk2 .* Ap2)   ;
    xk1   = xk1 + alpha .* pk1;
    xk2   = xk2 + alpha .* pk2;
    rk1p  = rk1 - alpha .* Ap1;
    rk2p  = rk2 - alpha .* Ap2;
    beta = ( rk1p.^2 + rk2p.^2 ) / ( rk1.^2 + rk2.^2 );
    pk1  = rk1p + beta * pk1;
    pk2  = rk2p + beta * pk2;
    
    rk1 = rk1p;
    rk2 = rk2p;
    
    % Second Step
    Ap1   = pk1 .* R1_dL + pk2 .* R1_dG;
    Ap2   = pk1 .* R2_dL + pk2 .* R2_dG;
    alpha = (rk1.^2 + rk2.^2)./(pk1 .* Ap1 + pk2 .* Ap2)   ;
    xk1   = xk1 + alpha .* pk1;
    xk2   = xk2 + alpha .* pk2;
    rk1   = rk1 - alpha .* Ap1;
    rk2   = rk2 - alpha .* Ap2;
%     beta = ( rk1p.^2 + rk2p.^2 ) / ( rk1.^2 + rk2.^2 );
%     pk1  = rk1p + beta * pk1;
%     pk2  = rk2p + beta * pk2;

    
    % Assign system solution
    dx    = [xk1,xk2];
    RNorm = abs(R1) + abs(R2);

end


