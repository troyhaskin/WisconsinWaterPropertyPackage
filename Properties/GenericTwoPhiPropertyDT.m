function Psi = GenericTwoPhiPropertyDT(delta,tau,OnePhiHandle,TwoPhiOption,PhaseCheck)
    
    if nargin < 3
        PhaseCheck = true;
    end
    
    Zeroer = delta*0;
    if (numel(tau) > numel(delta))
        Zeroer = tau*0;
    end
    
    Psi   = Zeroer*0                ; % Allocate
    
    if PhaseCheck && any(tau > 1)
        BelowCritical = (tau > 1);
        [Psat,dell,delg] = SaturationStateGivenTausat(tau(BelowCritical));
        TwoPhase = GetTwoPhaseMask(delta,dell,delg);

    else
        TwoPhase = false(size(Zeroer));
    end
    OnePhase = not(TwoPhase);
    
    
    if any(OnePhase)
        deltaOne      = SmartMask(delta,OnePhase);
        tauOne        = SmartMask(tau,OnePhase);
        Psi(OnePhase) = OnePhiHandle(deltaOne,tauOne,OnePhase);  
    end
    
    if any(TwoPhase)
        
        Choice = GetTwoPhiHandlingChoice(TwoPhiOption);
        
        switch(lower(Choice))
            case 'quality'
                delTwo  = SmartMask(delta,TwoPhase);
                dell    = SmartMask(dell ,TwoPhase);
                delg    = SmartMask(delg ,TwoPhase);
                tauTwo  = SmartMask(tau  ,TwoPhase);
                
                x    = QualityFromDensity(delTwo,dell,delg);
                Phil = OnePhiHandle(dell,tauTwo,TwoPhase);
                Phig = OnePhiHandle(delg,tauTwo,TwoPhase);
                
                Psi(TwoPhase) = Phil + x.*(Phig - Phil);
                
            case 'custom handle'
                Psi(TwoPhase) = TwoPhiOption(dell,delg,Psat,delta,tau,TwoPhase);
            otherwise
                
                error('Improper two-phase option/handle given.');
        end
    end
    
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

