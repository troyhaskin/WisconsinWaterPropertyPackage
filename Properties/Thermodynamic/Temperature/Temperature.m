function [T,varargout] = Temperature(rhoMix,iMix,Tguess,PhaseCheck)
    
    if (nargin < 4)
        PhaseCheck = true;
    end
    
    [rhoMix,SizeRho,iMix,SizeI] = Columnify(rhoMix,iMix)        ;
    [rhoMix, iMix]              = BalanceSizes(rhoMix,iMix)     ;
    
    if (nargin > 2) && all(not(isempty(Tguess))) && all(not(isnan(Tguess)))
        [rhoMix, Tguess] = BalanceSizes(rhoMix,Tguess)   ;
    end
    
    Size  = size(rhoMix);
    Sizer = ones(max(length(rhoMix),length(iMix)),1)    ;
    T     = Sizer * 0                                   ;
        
% 
% For a given mixture density, if the mixture internal energy
% is greater than the saturated internal energy, it is outside the vapor dome;
% therefore, the trouble of solving the two-phase mixture system can avoided.
%
    if PhaseCheck
        [~,Tsat,~,~] = SaturationStateGivenDensity(rhoMix) ; % Get Tsat from rhomix

        isat           = InternalEnergy(rhoMix,Tsat,false)   ; % Get isat
        BelowVaporDome = iMix < isat                         ;

        if any(BelowVaporDome)
            [rhoLt,rhoGt] = TriplePointDensities()                  ;
            Quality       = QualityFromDensity(rhoMix,rhoLt,rhoGt)  ;
            [itL,itG]     = TriplePointInternalEnergies()           ;
            itMix         = itL + Quality.*(itG - itL)              ;
            
            AboveTripleLine = iMix > itMix;
        else
            AboveTripleLine = true(Size); % If it is above the vapor dome, it is above the triple line
        end
        
    else
        BelowVaporDome  = false(Size)   ;
        AboveTripleLine = true(Size)    ;
    end

    % Phase masks
    OnePhase = not(BelowVaporDome) & AboveTripleLine;
    TwoPhase =     BelowVaporDome  & AboveTripleLine; % TwoPhase otherwise
    SolidGas = not(OnePhase) & not(TwoPhase);
    
    % begin: One-Phase Handling
    if any(OnePhase)
        
        % Pull One-phase values
        iOne        = SmartMask(iMix  ,OnePhase);
        rhoOne      = SmartMask(rhoMix,OnePhase);

        % Initial guess
        if ((nargin < 3) || isempty(Tguess)) && PhaseCheck
            TOne = SmartMask(1.1*Tsat,OnePhase);
            
        elseif (nargin < 3) || any(isempty(Tguess))
            TOne = 300*Sizer;
            
        else
            TOne = Tguess(OnePhase);
        end
        
        % Solve 
        T(OnePhase) = SolveOnePhase(rhoOne,iOne,TOne);
    end
    % end: One-Phase Handling
    
    
    % begin: Two-Phase Handling
    if any(TwoPhase)
        
        % Pull two-phase values
        iTwo        = SmartMask(iMix    ,TwoPhase);
        rhoTwo      = SmartMask(rhoMix  ,TwoPhase);
        
        % User-supplied initial guess
        if (nargin < 3) || isempty(Tguess)
            Ttwo = Tsat;
        else 
            Ttwo = Tguess(TwoPhase);
        end
        
        % Solve
        T(TwoPhase) = SolveTwoPhase(rhoTwo,iTwo,Ttwo,Tsat(TwoPhase));
    end
    % end: Two-Phase Handling


    % begin: Solid-Gas Phase Exception
    if any(SolidGas)
        
        warning('WWPP:Temperature:BelowTripleLine',...
                ['The state defined by the given density and internal energy fall ',...
                 'below the triple line which is an unsupported regime; therefore, ',...
                 'the temperature returned is the triple point temperature.']);
             
        T(SolidGas) = TriplePointTemperature();

    end
    % end: Solid-Gas Phase Exception

    T = RestoreShape(T,GreatestProduct(SizeRho,SizeI));
    
    if (nargout == 2)
        varargout{1} = any(TwoPhase);
    end
end

function T = SolveOnePhase(rho,i,Tguess)
    
    Tolerance = 1E-10;
    MaxIter   = 1E3     ;
    Updater   = @(T,Mask) Update(T,Mask,rho,i);
    T         = NewtonUpdater(Updater,Tguess,Tolerance,MaxIter);
    
    
    function [dT,Norm] = Update(T,Mask,rhoGiven,iGiven)
        iCalc     = InternalEnergy(rhoGiven(Mask),T,false)              ;
        Residual  = iGiven(Mask) - iCalc                                ;
        dResidual = -InternalEnergy_Temperature(rhoGiven(Mask),T,false) ;
        
        dT   = Residual ./ dResidual;
        Norm = abs(dT)        ;
    end
    
end

function T = SolveTwoPhase(rho,i,varargin)
    [~,T,~,~,~] = SaturationStateGivenMixedRhoI(rho,i); %,Tguess,Tsat);
end
