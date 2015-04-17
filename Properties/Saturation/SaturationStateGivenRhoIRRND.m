function varargout = SaturationStateGivenRhoIRRND(del,iND,varargin)
    
    [~,tauGuess,~,~] = SaturationStateGivenDeltaRRND(del);
    
    % Iteration Setup
    Guess     = tauGuess                ;
    IterMax   = 1E3                     ;
    Tolerance = 1E-12                   ;
    
    % Newton iterations
    [tau,Results] = NewtonClosure(@() MakeClosure(del,iND),Guess,Tolerance,IterMax);
    
    % Extract output values
    Pnd  = Results.Pnd;
    delL = Results.delL;
    delG = Results.delG;
    
    % Output
    Output    = {tau,Pnd,delL,delG} ;
    varargout = Output(1:nargout)   ;
    
end



function Closure = MakeClosure(del,iND)
    
    % Size of input
    N = length(del);
    
    % Initialize index filter
    Filter = (1:N)' ;
    
    % Initialize potential saturation quantities to return
    Zero = del * 0  ;
    
    % Initialize immutable arrays
    Pnd   = Zero;
    delL  = Zero;
    delG  = Zero;
    
    % Initialize values that will be reused; these are mutable.
    iL_      = Zero;
    iG_      = Zero;
    iLG_     = Zero;
    Pnd_     = Zero;
    delL_    = Zero;
    delG_    = Zero;
    x_       = Zero;
    del_     = del ;
    iND_     = iND ;
    OnePhase = true(N,1);
    TwoPhase = true(N,1);
    
    % Return the closure
    Closure.SetFilter    = @(x) SetFilter(x);
    Closure.GetResidual  = @(x) GetResidual(x);
    Closure.GetDResidual = @(x) GetResidualDerivative(x);
    Closure.Finalize     = @(x) Finalize (x);
    
    
    % ============================================================= %
    %       Setter for main Filter provided by NewtonClosure()      %
    % ============================================================= %
    function [] = SetFilter(Value)
        Filter = Value;
    end
    
    
    
    % ============================================================= %
    %        Getter for Residual requested by NewtonClosure()       %
    % ============================================================= %
    function R = GetResidual(tau)
        
        % Update the saturation state of the system and filter appropriately
        %         [Pnd_,delL_,delG_] = SaturationStateGivenTauRRND(tau);
        [Pnd,delL,delG] = AssignWithFilter(@()SaturationStateGivenTauRRND(tau),...
            Filter,Pnd,delL,delG);
        [Pnd_,delL_,delG_,del_,iND_] = FilterList(Filter,Pnd,delL,delG,del,iND);
        
        % Phase determination
        CannotSaturate = (   0  == Pnd_ )   ;     % Indicates it cannot saturate
        NotSaturated   = ( del_ < delG_ )   | ... % Super-heated gas
                         ( del_ > delL_ )   ;     % Sub-cooled liquid
        OnePhase       = CannotSaturate | NotSaturated  ;
        TwoPhase       = not(OnePhase)                  ;
        
        % Allocation
        R = tau;
        
        if any(OnePhase)
            R1  = InternalEnergyOneRND(del_(OnePhase),tau(OnePhase)) - iND_(OnePhase);
            
            % Assignment
            R (OnePhase) = R1   ;
        end
        
        if any(TwoPhase)
            
            [del2,tau2] = FilterList([TwoPhase;TwoPhase],[delL_;delG_],[tau;tau]);
            
            % Latent Internal Energy
            iLandG = InternalEnergyOneRND(del2,tau2);
            iL_     = iLandG(1:(end/2));
            iG_     = iLandG(end/2+1:end);
            iLG_    = iG_ - iL_;
            
            % Quality
            Masked = FilterList(TwoPhase,del_,delL_,delG_);
            x_ = QualityFromDensity(Masked{:});
            
            % Residual
            R2   = iL_ + x_.* iLG_ - iND(TwoPhase);
            
            % Assignment
            R (TwoPhase) = R2   ;
        end
        
        R = R./iND_;
    end
    
    
    
    % ============================================================= %
    %  Getter for Residual derivative requested by NewtonClosure()  %
    % ============================================================= %
    function dR = GetResidualDerivative(tau)
        
        % Filter
        [Pnd_,delL_,delG_,del_,iND_] = FilterList(Filter,Pnd,delL,delG,del,iND);
        
        % Phase determination
        CannotSaturate = (   0  == Pnd_ )   ;     % Indicates it cannot saturate
        NotSaturated   = ( del_ < delG_ )   | ... % Super-heated gas
            ( del_ > delL_ )   ;     % Sub-cooled liquid
        OnePhase       = CannotSaturate | NotSaturated  ;
        TwoPhase       = not(OnePhase)                  ;
        
        % Allocation
        dR = Pnd_;
        
        % One phase residual
        if any(OnePhase)
            dR1 = InternalEnergyOneRND_tau(del_(OnePhase),tau(OnePhase));
            
            % Assignment
            dR (OnePhase) = dR1   ;
        end
        
        
        % Two phase residual
        if any(TwoPhase)
            
            % Filter
            [tau2,Pnd2,delL2,delG2,del2] = FilterList(TwoPhase,tau,Pnd_,delL_,delG_,del_);
            
            delLandG = [delL2;delG2];
            tauLandG = [tau2;tau2];
            FilterL  = 1:length(delL2);
            FilterG  = length(delL2)+1 : 2*length(delL2);
            
            % Helmholtz calls
            PhiRLG_d  = HelmholtzResidual_d (delLandG,tauLandG);
            PhiRLG_dd = HelmholtzResidual_dd (delLandG,tauLandG);
            PhiRLG_dt = HelmholtzResidual_dt (delLandG,tauLandG);
            diLandG_d = InternalEnergyOneRND_delta(delLandG,tauLandG);
            diLandG_t = InternalEnergyOneRND_tau(delLandG,tauLandG);
            
            % Derivative of iL w.r.t. tau
            Pnd_tau  = ClausiusClapeyronRRND(Pnd2,tau2,delL2,delG2,iL_,iG_);
            PhiR_d   = PhiRLG_d (FilterL);
            PhiR_dd  = PhiRLG_dd(FilterL);
            PhiR_dt  = PhiRLG_dt(FilterL);
            delL_tau = (delL2 + tau2.^2 .* Pnd_tau + delL2.^2 .* (PhiR_d - tau2 .* PhiR_dt)) ./ ...
                (tau2 .* (1 + 2 * delL2 .* PhiR_d + delL2.^2 .* PhiR_dd ));
            diL = diLandG_d(FilterL) .* delL_tau + diLandG_t(FilterL);
            
            % Derivative of iG w.r.t. tau
            PhiR_d   = PhiRLG_d (FilterG);
            PhiR_dd  = PhiRLG_dd(FilterG);
            PhiR_dt  = PhiRLG_dt(FilterG);
            delG_tau = (delG2 + tau2.^2 .* Pnd_tau + delG2.^2 .* (PhiR_d - tau2 .* PhiR_dt)) ./ ...
                (tau2 .* (1 + 2 * delG2 .* PhiR_d + delG2.^2 .* PhiR_dd ));
            diG = diLandG_d(FilterG) .* delL_tau + diLandG_t(FilterG);
            
            % Derivative of iLG w.r.t. tau
            diLG = diG - diL;
            
            % Derivative of x w.r.t. tau
            dx = ((delL2-del2).*delL2.*delG_tau + (del2-delG2).*delG2.*delL_tau)./(del2.*(delG2-delL2).^2);
            
            % Quality
            x_ = QualityFromDensity(del2,delL2,delG2);
            
            % Derivative of residual w.r.t. tau
            dR2 = diL + x_ .* diLG + (iG_ - iL_) .* dx;
            
            % Assignment
            dR (TwoPhase) = dR2   ;
        end
        
        dR = dR./iND_;
    end
    
    
    
    % ============================================================= %
    %       Getter for finalization requested by NewtonClosure()    %
    % ============================================================= %
    function Results = Finalize(tau)
        
        [Pnd,delL,delG] = SaturationStateGivenTauRRND(tau);
        
        Results.tau  = tau  ;
        Results.delL = delL ;
        Results.delG = delG ;
        Results.Pnd  = Pnd  ;
    end
    
    
    
end

