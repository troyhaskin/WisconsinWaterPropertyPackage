function [Psat,tauSat,delL,delG] = SaturationStateGivenDeltaRRND(delKn,UniqueMask)
    
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
    
    Solution = NewtonUpdater(@(Unknowns,Mask) UpdateSystem(Unknowns,Mask,p,mu,delL),tauGuess,1E-14,100,true);
end

function [dUnknowns,RNorm] = UpdateSystem(Unknowns,Mask,delGiven,N)
    
    % Pull unknowns
    dels = Unknowns(:,1);
    taus = Unknowns(:,2);

    % Generate masks for the Unknowns
    MaskL = (Mask <= N) ;
    MaskG = not(MaskL)  ;
    
    % delG given for saturation solution
    delLL = dels    (MaskL);
    tauL  = taus    (MaskL);
    delGL = delGiven(MaskL);

    % delL given for saturation solution
    delGG = dels    (MaskG);
    tauG  = taus    (MaskG);
    delLG = delGiven(MaskG);
    
    % Combined vectors for liquid/gas-given agnostic (LGA) evaluations (aides in vectorization)
    delL = [delLL ; delLG ];
    delG = [delGL ; delGG ];
    tau  = [tauL  ; tauG  ];


    % Helmholtz energy function handles (d = del, t = tau; shortened to not 
    % conflict with above variables).
    PhiR    = @(d,t) HelmholtzResidual   (d,t);
    PhiR_d  = @(d,t) HelmholtzResidual_d (d,t);
    PhiR_t  = @(d,t) HelmholtzResidual_t (d,t);
    PhiR_dd = @(d,t) HelmholtzResidual_dd(d,t);
    PhiR_dt = @(d,t) HelmholtzResidual_dt(d,t);


    % ============================================================================ %
    %                      Calculate Dimensionless Pressure                        %
    % ============================================================================ %
    % Pnd is the dimensionless pressure eliminated from
    % the system to shrink the solution space.  Here are the functions and
    % derivatives that are needed in the residuals below.

    % Pnd = Pnd1 * Pnd2 to more easily maintain/check the handles.
    Pnd1 = PhiR(delL,tau) - PhiR(delG,tau) + log(delL./delG); % LGA
    Pnd2 = (delL.*delG)./(delL - delG)./tau                 ; % LGA

    % tau derivatives  (LGA)
    Pnd1_t = PhiR_t(delL,tau) - PhiR_t(delG,tau);
    Pnd2_t = -Pnd2./tau                         ;

    % del derivatives while solving for delL, tauL with delG given
    if any(MaskL)
        Pnd1_dL = 1./delLL + PhiR_d(delLL,tauL)  ;
        f       = delGL./(delLL - delGL)./tauL   ;
        Pnd2_dL = (1 - delLL./(delLL-delGL)).*f  ;
    else
        Pnd1_dL = [];
        Pnd2_dL = [];
    end
    
    % del derivatives while solving for delG, tauG with delL given
    if any(MaskL)
        Pnd1_dG = -1./delGG + PhiR_d(delGG,tauG)  ;
        f        = delLG./(delLG - delGG)./tauG   ;
        Pnd2_dG = (1 + delLG./(delLG-delGG)).*f  ;
    else
        Pnd1_dG = [];
        Pnd2_dG = [];
    end
    
    % Create combined Pnd*_d vectors
    Pnd1_d = [Pnd1_dL ; Pnd1_dG];
    Pnd2_d = [Pnd2_dL ; Pnd2_dG];
    
    % Pnd calculated
    Pnd    = Pnd1   .* Pnd2                    ;
    Pnd_d  = Pnd1_d .* Pnd2 + Pnd2_d .* Pnd1   ;
    Pnd_t  = Pnd1_t .* Pnd2 + Pnd2_t .* Pnd1   ;

    % Residuals
    R1   = Pnd  - delL .* (1 + delL.*PhiR_d(delL,tau));
    R2   = Pnd  - delG .* (1 + delG.*PhiR_d(delG,tau));
    
    % Jacobian values
    R1_d = Pnd_d - (1 + 2*delL.*PhiR_d(delL,tau) + delL.^2.*PhiR_dd(delL,tau))  ;
    R1_t = Pnd_t - delL.^2 .* PhiR_dt(delL,tau)                                 ;
    R2_d = Pnd_d                                                                ;
    R2_t = Pnd_t - delG.^2 .* PhiR_dt(delG,tau)                                 ;
    DetR = R1_d .* R2_t - R1_t.*R2_d                                            ;


    % Final Newton Directions
    ddel  =  (R1 .* R2_t - R2 .* R1_t) ./ DetR ;
    dtau  =  (R2 .* R1_d - R1 .* R2_d) ./ DetR ;

    % Final outputs
    dUnknowns = [ddel,dtau]         ;
    RNorm     = abs(R1) + abs(R2)   ; % L_1 norm

end