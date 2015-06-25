function Psi = GenericTwoPhiPropertyR(delta,tau,OnePhiHandle,TwoPhiOption,PhaseCheck,twoPhiState)
    
    if nargin < 5
        PhaseCheck = true;
    end
    
    %     Allocation and Setup
    % ----------------------------
    Nstates = max(length(delta),length(tau));
    Nprops  = size(OnePhiHandle(NaN,NaN,true(Nstates,1)),2);
    Psi     = zeros(Nstates,Nprops);


    %  Phase Determination
    % -----------------------------------------------
    if isempty(twoPhiState) || all(not(twoPhiState.isTwoPhi))
        % Perform the phase boundary check
        if PhaseCheck
            [Pnd,delL,delG] = SaturationStateGivenTauRRND(tau);
            isTwoPhi = GetTwoPhaseMask(delta,delL,delG);
        else
            isTwoPhi = false(Nstates,1);
        end
    else
        delL     = twoPhiState.delL     ;
        delG     = twoPhiState.delG     ;
        Pnd      = twoPhiState.Pnd      ;
        isTwoPhi = twoPhiState.isTwoPhi ;
    end
    isOnePhi = not(isTwoPhi);
    %   end: Phase Determination
    % -----------------------------------------------
    
    
    % -----------------------------------------------
    % begin: One Phase Handling
    if any(isOnePhi)
        del1 = SmartMask(delta,isOnePhi);
        tau1 = SmartMask(tau  ,isOnePhi);
        
        Psi(isOnePhi,:) = OnePhiHandle(del1,tau1,isOnePhi);
    end
    % end: One Phase Handling
    
    
    % begin: Two Phase Handling
    if any(isTwoPhi)
        
        % Grab the two-phase values only
        del  = SmartMask(delta,isTwoPhi);
        delL = SmartMask(delL ,isTwoPhi);
        delG = SmartMask(delG ,isTwoPhi);
        Pnd  = SmartMask(Pnd  ,isTwoPhi);
        tau  = SmartMask(tau  ,isTwoPhi);
        
        % Determine the two-phase methodology passed in
        Choice  = GetTwoPhiHandlingChoice(TwoPhiOption);
        
        % begin: Generation of Two-phase Properties
        switch(lower(Choice))
            
            case('quality') % Quality-weighted
                x    = QualityFromDensity(del,delL,delG) ; % Quality
                PsiL = OnePhiHandle(delL,tau,isTwoPhi)    ; % Sat. Liquid value
                PsiG = OnePhiHandle(delG,tau,isTwoPhi)    ; % Sat. Gas value
                
                LatentPsi       = PsiG - PsiL;
                Psi(isTwoPhi,:) = PsiL + bsxfun(@times,x,LatentPsi); % Quality-weighted value
                
                
            case({'void','void fracion','homogeneous void fraction'}) % Void-weighted
                
                alpha = HomogeneousVoidFraction(del,delL,delG)   ; % Homogeneous Void Fraction
                PsiL  = OnePhiHandle(delL,tau,isTwoPhi)           ; % Sat. Liquid value
                PsiG  = OnePhiHandle(delG,tau,isTwoPhi)           ; % Sat. Gas value
                
                LatentPsi       = PsiG - PsiL;
                Psi(isTwoPhi,:) = Psil + bsxfun(@times,alpha,LatentPsi); % Void-weighted value
                
                
            case('customhandle') % Custom weight function
                Psi(isTwoPhi,:) = TwoPhiOption(Pnd,delL,delG,del,tau,isTwoPhi);
                
            case('undefined')
                warning('Thermodynamics:GenericTwoPhiPropertyR:UndefinedTwoPhaseBehavior',...
                        'Two phase behavior for the property is undefined.');

            otherwise
                error('Thermodynamics:GenericTwoPhiPropertyR:UnspecifiedTwoPhaseHandling',...
                      'Improper two-phase option ''%s'' given.',Choice);

        end
        % end: Generation of Two-phase Properties
    end
    % end: Two Phase Handling

end



function Choice = GetTwoPhiHandlingChoice(TwoPhiOption)
    if ischar(TwoPhiOption)
        Choice = TwoPhiOption;
        
    elseif isa(TwoPhiOption,'function_handle')
        Choice = 'CustomHandle';
        
    else
        Choice = '';
        
    end
end


function TwoPhase = GetTwoPhaseMask(rho,rhol,rhog)
    SizesMatch = size(rho) == size(rhol);
    
    if all(SizesMatch) || isscalar(rhol) % do straight vector comparison
        TwoPhase = (rho < rhol) & (rho > rhog);
        
    elseif any(SizesMatch) % perform singleton expansion with a transpose
        BelowLiquidDensity = bsxfun(@lt,rho,rhol');
        AboveGasDensity    = bsxfun(@gt,rho,rhog');
        
        TwoPhase = BelowLiquidDensity & AboveGasDensity;
        
        if not(isvector(TwoPhase))
            TwoPhase = any(TwoPhase);
        end
        
    else % perform singleton expansion without a transpose
        BelowLiquidDensity = any(bsxfun(@lt,rho,rhol));
        AboveGasDensity    = any(bsxfun(@gt,rho,rhog));
        
        TwoPhase = BelowLiquidDensity & AboveGasDensity;
    end
end

