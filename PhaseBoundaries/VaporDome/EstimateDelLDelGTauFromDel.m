function [delL,delG,tau] = EstimateDelLDelGTauFromDel(del)
    
    % Allocate
    delL = del * 0;
    delG = del * 0;
    tau  = del * 0;
    
    % Phase logicals
    IsDelC   = del == 1                 ;   % critical state
    IsDelL   = del >  1                 ;   % liquid side in stable Newton iteration regime
    IsDelG   = del <  1                 ;   % gas    side in stable Newton iteration regime
    IterateL = IsDelL & (del > 1.0005)  ;   % Stable Newton Iteration regime
    IterateG = IsDelG & (del < 0.9984)  ;   % Stable Newton Iteration regime
    
    % Push known values
    delL(IsDelL) = del(IsDelL);
    delG(IsDelG) = del(IsDelG);
    
    % Catch critical states and auto-assign (also avoids need for iteration)
	delL(IsDelC) = 1    ;
    delG(IsDelC) = 1    ;
    tau (IsDelC) = 1    ;
    
    % Iteration parameters
    Tolerance = DefaultAbsoluteIterationTolerance() ;
    IterMax   = DefaultMaximumIterationCount()      ;
    
    % Solve with liquid density
    if any(IsDelL)
        tauL = EstimateTauFromDelL(del(IsDelL));                 % Get initial guess

        if any(IterateL)
            Updater = @(t,Mask) GetTauFromDelL(t,Mask,del(IterateL)); % Update handle
            tauL(IterateL(IsDelL)) = NewtonUpdater(Updater,tauL(IterateL(IsDelL)),Tolerance,IterMax);   % Solve
        end

        tau (IsDelL) = tauL    ; % Assign tau value
        delG(IsDelL) = EstimateDelGFromTau(tauL) ; % Assign delG value
    end
    
    % Solve with gas density
    if any(IsDelG)
        tauG = EstimateTauFromDelG(del(IsDelG));                 % Get initial guess
        
        if any(IterateG)
            Updater = @(t,Mask) GetTauFromDelG(t,Mask,del(IterateG)); % Update handle
            tauG(IterateG(IsDelG))    = NewtonUpdater(Updater,tauG(IterateG(IsDelG)),Tolerance,IterMax);   % Solve
        end

        tau (IsDelG) = tauG    ; % Assign tau value
        delL(IsDelG) = EstimateDelLFromTau(tauG) ; % Assign delL value
    end

end



function [dx,Norm] = GetTauFromDelL(tau,Mask,delL)
    
    Residual  = EstimateDelLFromTau(tau) - delL(Mask);
    dResidual = EstimateDelLFromTau_tau(tau);

    dx      = Residual ./ dResidual ;
    Norm = abs(dx);
end

function [dx,Norm] = GetTauFromDelG(tau,Mask,delG)
    
    Residual  = EstimateDelGFromTau(tau) - delG(Mask);
    dResidual = EstimateDelGFromTau_tau(tau);
    
    dx   = Residual ./ dResidual;
    Norm = abs(Residual)        ;
    
end
% 
% 
% function tau = EstimateTauFromNearCriticalDelL(delL)
%     p = [   -1.2427888617613720E-13 ,...
%             +1.6349800613368345E-11 ,...
%             -1.6826603201245960E-09 ,...
%             +1.7914607206586920E-07 ,...
%             +1.5044848361890351E-06 ,...
%             +4.0702911932209303E-06 ,...
%             +1.0000036400293815E+00 ];
%     mu = [+1.0309140731308570E+00,+1.1619681188514367E-02];
%     
%     tau = HornersMethod(delL,p,mu);
% end
% 
% function tau = EstimateTauFromNearCriticalDelG(delG)
%     p = [   -5.8888931624887446E-13 ,...
%             -3.3590363359348875E-11 ,...
%             -5.1209232832641621E-10 ,...
%             -1.8576024219530422E-07 ,...
%             +1.5151724268411565E-06 ,...
%             -4.0702210622798952E-06 ,...
%             +1.0000036325648092E+00 ];
%     mu = [+9.6863921202293113E-01,+1.1757101679711156E-02];
%     
%     tau = HornersMethod(delG,p,mu);
% end
% 
% 
