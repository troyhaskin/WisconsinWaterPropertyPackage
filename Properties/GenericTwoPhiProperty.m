function Psi = GenericTwoPhiProperty(rho,T,OnePhiHandle,TwoPhiOption,PhaseCheck)
       
    if nargin < 5
        PhaseCheck = true;
    end
    
    % -----------------------------------------------
    % begin: Allocation and Setup
    Nstates = max(length(rho),length(T));
    Nprops  = size(OnePhiHandle(NaN,NaN,true(Nstates,1)),2);
    Psi     = zeros(Nstates,Nprops);
    
    [~,rhoc,Tc] = Nondimensionalizers();
    [rho,T]     = BalanceSizes(rho,T);
    % end:  Allocation and Setup
    
    
    % -----------------------------------------------
    % begin: Reduce the Quantities
    
    tau = Tc ./ T   ;
    del = rho / rhoc;
    
    % end: Reduce the Quantities

    
    % -----------------------------------------------
    %   begin: Phase Determination
    if PhaseCheck
        [Psat,delL,delG] = SaturationStateGivenTausat(tau)  ;
        TwoPhase         = (del < delL) & (del > delG)      ;
    else
        TwoPhase = false(Nstates,1);
    end
    OnePhase = not(TwoPhase);
    %   end: Phase Determination
    % -----------------------------------------------

    
    % -----------------------------------------------
    % begin: One Phase Handling
    if any(OnePhase)
        del1 = SmartMask(del,OnePhase);
        tau1 = SmartMask(tau,OnePhase);
        
        Psi(OnePhase,:) = OnePhiHandle(del1,tau1,OnePhase);  
    end
    % end: One Phase Handling
    
    
    % begin: Two Phase Handling
    if any(TwoPhase)
        
        % Grab the two-phase values only
        del2    = SmartMask(del ,TwoPhase)  ;
        delL2   = SmartMask(delL,TwoPhase)  ;
        delG2   = SmartMask(delG,TwoPhase)  ;
        Psat2   = SmartMask(Psat,TwoPhase)  ;
        tau2    = SmartMask(tau ,TwoPhase)  ;
        
        % Determine the two-phase methodology passed in
        Choice  = GetTwoPhiHandlingChoice(TwoPhiOption);
        
        % begin: Generation of Two-phase Properties
        switch(lower(Choice))
            
            case('quality') % Quality-weighted
                x    = QualityFromDensity(del2,delL2,delG2) ; % Quality
                Psil = OnePhiHandle(delL2,tau2,TwoPhase)    ; % Sat. Liquid value
                Psig = OnePhiHandle(delG2,tau2,TwoPhase)    ; % Sat. Gas value
                
                LatentPsi       = Psig - Psil;
                Psi(TwoPhase,:) = Psil + bsxfun(@times,x,LatentPsi); % Quality-weighted value


            case({'void','void fracion','homogeneous void fraction'}) % Void-weighted
                
                alpha = HomogeneousVoidFraction(del2,delL2,delG2)   ; % Homogeneous Void Fraction
                Psil  = OnePhiHandle(delL2,tau2,TwoPhase)           ; % Sat. Liquid value
                Psig  = OnePhiHandle(delG2,tau2,TwoPhase)           ; % Sat. Gas value
                
                LatentPsi       = Psig - Psil;
                Psi(TwoPhase,:) = Psil + bsxfun(@times,alpha,LatentPsi); % Void-weighted value


            case('custom handle') % Custom weight function
                Psi(TwoPhase,:) = TwoPhiOption(Psat2,delL2,delG2,del2,tau2,TwoPhase);
                
            otherwise
                
                error('Thermodynamics:GenericTwoPhiProperty:BadTwoPhaseOption',...
                      'Improper two-phase option/handle given.');
        end
        % end: Generation of Two-phase Properties
    end
    % end: Two Phase Handling
    
end



function Choice = GetTwoPhiHandlingChoice(TwoPhiOption)
    if ischar(TwoPhiOption)
        Choice = TwoPhiOption;

    elseif isa(TwoPhiOption,'function_handle')
        Choice = 'Custom Handle';

    else
        Choice = '';

    end
end

