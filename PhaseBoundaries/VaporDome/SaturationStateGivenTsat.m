function [Psat,rhol,rhog] = SaturationStateGivenTsat(Tsat,rhol0,rhog0)
    
    [Tsat,SizeTsat] = Columnify(Tsat);
    [Tsat,~,UniqueMask] = unique(Tsat);
    
    Tc   = CriticalTemperature()    ; %[K]
    rhoc = CriticalDensity()        ; %[kg/m^3]
    tau  =  Tc ./ Tsat              ; %[-] - inverse reduced temperature
    
    
    % Reduced liquid density guess value
    if (nargin < 2)
        delL  = EstimateDelLFromTau(tau);
    else
        delL = rhol0 / rhoc;
    end
	rhol  = delL * rhoc;
    
    % Reduced gas density guess value
    if (nargin < 2)
        delG  =  EstimateDelGFromTau(tau);
    else
        delG = rhog0 / rhoc;
    end
	rhog  = delG * rhoc;
       

    % Check Tsat for how close to critical it is
    NearTc  = tau < DoNoIterationValueReducedTemperature()  ; 
    AboveTc = tau <   1                                     ; % Check if Tsat is above the critical point

    % Iteration setup
    Tolerance = 1E-12                                       ; % Abolsute iteration tolerance
    IterMax   = 1E3                                         ; % Maximum iteration count
    Calculate = not(NearTc) & not(AboveTc)                  ; % Mask for temps. not close to the critical point
    Guess     = [delL(Calculate),delG(Calculate)]           ; % Starting values for the iteration
    UpdateFun = @(x,Mask) Updater(x,Mask,tau(Calculate))    ; % Function used for updating the solution
    ObjectiveFun = @(x,Mask) Objective(x,Mask,tau(Calculate));
    
    % Solve the system
    if any(Calculate)
        %xSol = NewtonUpdater(UpdateFun,Guess,Tolerance,IterMax);
        xSol = BackTrackingNewtonUpdater(UpdateFun,ObjectiveFun,Guess,Tolerance,IterMax);
        
        % Update the iterated values
        rhol(Calculate) = xSol(:,1) * rhoc;
        rhog(Calculate) = xSol(:,2) * rhoc;
    end
    
    % These values above for T above the vapor dome ensure that the fluid will
    % be flagged as single phase.
    rhol(AboveTc) = 0.0;
    rhog(AboveTc) = 1E3;
    
    % Expand vectors back to non-unique lengths
    Tsat = Tsat(UniqueMask)         ;
    rhol = rhol(UniqueMask)         ;
    rhog = rhog(UniqueMask)         ;
    
    % Calculate the saturation pressure
    Psat = Pressure(rhol,Tsat,false);
    
    % Reshape to input shape
    Psat = RestoreShape(Psat,SizeTsat);
    rhol = RestoreShape(rhol,SizeTsat);
    rhog = RestoreShape(rhog,SizeTsat);
    
end


function [dx,RNorm] = Updater(x,Mask,tau0)
    
    tau  = tau0(Mask);
    delL = x(:,1);
    delG = x(:,2);
    
    PhiR    = @(delta) HelmholtzResidual   (delta,tau);
    PhiR_d  = @(delta) HelmholtzResidual_d (delta,tau);
    PhiR_dd = @(delta) HelmholtzResidual_dd(delta,tau);
    
    Psig1    = PhiR(delL) - PhiR(delG) + log(delL./delG)    ;
    Psig1_dl =  PhiR_d(delL) + 1./delL                      ;
    Psig1_dg =-(PhiR_d(delG) + 1./delG)                     ;
    
    PdelL = 1 + delL.*PhiR_d(delL);
    PdelG = 1 + delG.*PhiR_d(delG);
    
    PdelL_d = PhiR_d(delL) + delL.*PhiR_dd(delL);
    PdelG_d = PhiR_d(delG) + delG.*PhiR_dd(delG);
    
    R1 = delG.*Psig1 - PdelL.*(delL-delG);
    R2 = delL.*Psig1 - PdelG.*(delL-delG);
    
    R1_dl = delG.*Psig1_dl + (delG-delL).*PdelL_d - PdelL;
    R1_dg = Psig1 + PdelL + delG.*Psig1_dg;
    
    R2_dl = Psig1 - PdelG + delL.*Psig1_dl;
    R2_dg = delL.*Psig1_dg + (delG-delL).*PdelG_d + PdelG;
    
    DetJ = R1_dl .* R2_dg - R1_dg .*R2_dl;
    
    iJ11 =  R2_dg ./ DetJ;
    iJ12 = -R1_dg ./ DetJ;
    iJ21 = -R2_dl ./ DetJ;
    iJ22 =  R1_dl ./ DetJ;
    
    ddelL = iJ11 .* R1 + iJ12 .* R2;
    ddelG = iJ21 .* R1 + iJ22 .* R2;
    
    dx    = [ddelL,ddelG];
    RNorm = abs(R1) + abs(R2);
    
end

function RNorm = Objective(x,Mask,tau0)
    
    if any(x(:) <= 0)
        RNorm = NaN;
        return
    end
    
    tau  = tau0(Mask);
    delL = x(:,1);
    delG = x(:,2);

    PhiR    = @(delta) HelmholtzResidual   (delta,tau);
    PhiR_d  = @(delta) HelmholtzResidual_d (delta,tau);
    
    Psig1    = PhiR(delL) - PhiR(delG) + log(delL./delG)    ;
    
    PdelL = 1 + delL.*PhiR_d(delL);
    PdelG = 1 + delG.*PhiR_d(delG);
    
    R1 = delG.*Psig1 - PdelL.*(delL-delG);
    R2 = delL.*Psig1 - PdelG.*(delL-delG);
    
    RNorm = abs(R1) + abs(R2);
end