function [Helm,Helm_d] = HelmholtzResidualCombo__d(delta,tau)
    
    [c,d,t,n,alpha,beta,gamma,epsilon,A,B,C,D,a,b] = Coefficients_HelmholtzResidual();

    
    deltaMod	= (delta - 1).^2 + eps(delta)    ;
    betaInv     = 1 ./ beta         ;
    Helm      = 0*delta               ;
    SumErr    = 0*delta               ; % Summation correction (Kahan summation)
    
    for k = 1:7
        %
        % Calculate zeroth derivative Helmholtz free energy component
        Part    = n(k) * delta.^(d(k)) .* tau.^(t(k))   ;
        [Helm,SumErr] = KahanSum(Helm,Part,SumErr)      ;
        %
        % Calculate first derivative Helmholtz free energy component
        Part    = n(k) * d(k) * delta.^(d(k)-1) .* tau.^(t(k))  ;
        [Helm_d,SumErr] = KahanSum(Helm_d,Part,SumErr)    ;
    end
    
    for k = 8:51
        %
        % Calculate zeroth derivative Helmholtz free energy component
        Part    = n(k) * delta.^(d(k)) .* tau.^(t(k)) .* exp(-delta.^(c(k)))    ;
        [Helm,SumErr] = KahanSum(Helm,Part,SumErr);
        %
        % Calculate first derivative Helmholtz free energy component
        Part    = d(k) - c(k) * delta.^(c(k))               ;
        Part    = delta.^(d(k)-1) .* tau.^(t(k)).*(Part)    ;
        Part    = n(k) * exp(-delta.^(c(k))) .* Part        ;
        [Helm_d,SumErr] = KahanSum(Helm_d,Part,SumErr);
    end
    
    for k = 52:54
        m       = k - 51                                                        ;
        Arg     = -alpha(m)*(delta-epsilon(m)).^2 - beta(m)*(tau-gamma(m)).^2   ;
        %
        % Calculate zeroth derivative Helmholtz free energy component
        Part    = n(k) * delta.^(d(k)) .* tau.^(t(k)) .* exp(Arg)   ;
        [Helm,SumErr] = KahanSum(Helm,Part,SumErr)                  ;
        %
        % Calculate first derivative Helmholtz free energy component
        Part    = n(k) * delta.^(d(k)) .* tau.^(t(k)) .* exp(Arg) .*(d(k)./delta - 2*alpha(m).*(delta-epsilon(m)));
        [Helm_d,SumErr] = KahanSum(Helm_d,Part,SumErr);
    end
    
    for k = 55:56
        m        = k - 51	;
        p        = k - 54	;

        Theta       = GetTheta     (deltaMod,tau,A(p),betaInv(m))                             ;
        Delta       = GetDelta     (deltaMod,Theta,B(p),a(p))                                 ;
        Psi         = GetPsi       (deltaMod,tau,C(p),D(p))                                   ;
        Psi_d       = GetPsi_d     (delta,Psi,C(p))                                           ;
        Deltabi_d   = GetDeltabi_d (delta,deltaMod,Delta,Theta,A(p),B(p),a(p),b(p),betaInv(m));
        %
        % Calculate zeroth derivative Helmholtz free energy component
        Part    = n(k) * Delta.^(b(p)) .* delta .* Psi  ;
        [Helm,SumErr] = KahanSum(Helm,Part,SumErr);
        %
        % Calculate first derivative Helmholtz free energy component
        Part        =  Delta.^(b(p)) .* (Psi + delta .* Psi_d)  ;
        Part        = n(k)*(Part + Deltabi_d .* delta .* Psi)   ;
        [Helm_d,SumErr] = KahanSum(Helm_d,Part,SumErr)    ;
    end
    
end
