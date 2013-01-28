function T = Temperature(rhoMix,iMix,Tguess)
    

% % Load scaling constants
%     Tc    = CriticalTemperature()           ; %[K]
%     rhoc  = CriticalTemperature()           ; %[kg/m^3]
    
% Allocate the solution vector
    switch(numel(rhoMix) > numel(iMix))
        case true
            T = zeros(size(rhoMix));
        case false
            T = zeros(size( iMix ));
    end
        
% 
% For a given mixture density, if the mixture internal energy
% is greater than the saturated internal energy, it is outside the vapor dome;
% therefore, the trouble of solving the two-phase mixture system can avoided.
%
    [~,Tsat,~,~] = SaturationStateGivenDensity(rhoMix) ; % Get Tsat from rhomix
    isat         = InternalEnergy(rhoMix,Tsat,false)   ; % Get isat
    OnePhase     = iMix > isat                         ; % OnePhase if imix is above isat
    TwoPhase     = not(OnePhase)                       ; % TwoPhase otherwise
    
    
    if any(OnePhase)
        iOne        = SmartMask(iMix    ,OnePhase);
        rhoOne      = SmartMask(rhoMix  ,OnePhase);

        if (nargin < 3)
            TOne = SmartMask(1.1*Tsat,OnePhase);
        else
            TOne = Tguess(OnePhase);
        end
        
        T(OnePhase) = SolveOnePhase(rhoOne,iOne,TOne);
    end
    
    if any(TwoPhase)
        iTwo        = SmartMask(iMix    ,TwoPhase);
        rhoTwo      = SmartMask(rhoMix  ,TwoPhase);
        
        if (nargin < 3)
            Tfreeze = 273.15;
            TTwo    = SmartMask((Tfreeze + Tsat)/2,TwoPhase); % In-the-middle is as good a guess as any
        else
            TTwo = Tguess(TwoPhase);
        end
        
        T(TwoPhase) = SolveTwoPhase(rhoTwo,iTwo,TTwo);
    end

end

function T = SolveOnePhase(rho,i,Tguess)
	
    rhoc  = CriticalDensity()               ; %[kg/m^3]
    Tc    = CriticalTemperature()           ; %[K]
    ic    = InternalEnergy(rhoc,Tc,false)   ; %[J/kg]
    
	iota  = @(T,Mask) InternalEnergy(rho(Mask),T)/ic;
	diota = @(T,Mask) InternalEnergy_Temperature(rho(Mask),T)/ic;
                
    R  = @(T,Mask) i(Mask)/ic - iota (T,Mask) ;
    dR = @(T,Mask)            - diota(T,Mask) ;
    
    Tolerance = 1E-12   ;
    MaxIter   = 1E3     ;
    T         = NewtonsMethodOneVariable(R,dR,Tguess,Tolerance,MaxIter);
end

function T = SolveTwoPhase(rho,i,Tguess)
    [~,T,~,~,~] = SaturationStateGivenRhoImix(rho,i,Tguess);
end
