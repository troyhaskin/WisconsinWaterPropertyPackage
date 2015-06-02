function [Helm,Helm_d] = HelmholtzResidualCombo__d(delta,tau)
    
    [c,d,t,n,alpha,beta,gamma,epsilon,A,B,C,D,a,b] = Coefficients_HelmholtzResidual();

    
    betaInv       = 1 ./ beta   ;
    N             = numel(delta);
    Helm(N,1)     = 0           ;
    Helm_d(N,1)   = 0           ;
%     SumErr(N,1)   = 0           ; % Summation correction (Kahan summation)
%     SumErr_d(N,1) = 0           ; % Summation correction (Kahan summation)
    
    for k = 1:7
        
        % Shared portion
        Part = n(k) * delta.^(d(k)-1) .* tau.^(t(k));
        
        
        % Calculate zeroth derivative Helmholtz free energy component
        %         [Helm,SumErr] = KahanSum(Helm, delta.^2 .* Part ,SumErr)      ;
        Helm = Helm + delta .* Part;
        
        
        % Calculate first derivative Helmholtz free energy component
        %         [Helm_d,SumErr_d] = KahanSum(Helm_d, d(k) * delta .* Part ,SumErr_d)      ;
        Helm_d = Helm_d + d(k) * Part;

    end
    
    for k = 8:51
        % Helper variables
%         bsxfun(@power,delta,c(k:k+1));
        delta2c = delta.^(c(k));
        
        % Calculate zeroth derivative Helmholtz free energy component
        Part    = n(k) * delta.^(d(k)) .* tau.^(t(k)) .* exp(-delta2c)  ;
        %         [Helm,SumErr] = KahanSum(Helm,Part,SumErr)                      ;
        Helm = Helm + Part;
        
        
        % Calculate first derivative Helmholtz free energy component
        Part_d = Part .* (d(k) - c(k) * delta2c) ./ delta   ;
        %         [Helm_d,SumErr_d] = KahanSum(Helm_d,Part_d,SumErr_d);
        Helm_d = Helm_d + Part_d;

    end
    
    for k = 52:54
        m       = k - 51                                                        ;
        Arg     = -alpha(m)*(delta-epsilon(m)).^2 - beta(m)*(tau-gamma(m)).^2   ;
        
        
        % Calculate zeroth derivative Helmholtz free energy component
        %         [Helm,SumErr] = KahanSum(Helm,Part,SumErr)                  ;
        Helm = Helm + n(k) * delta.^(d(k)) .* tau.^(t(k)) .* exp(Arg);
        
        
        % Calculate first derivative Helmholtz free energy component
        %         [Helm_d,SumErr_d] = KahanSum(Helm_d,Part_d,SumErr_d)            ;
        Helm_d = Helm_d + Part .* (d(k)./delta - 2*alpha(m)*(delta-epsilon(m)));

    end


    notCrit  = (delta ~= 1)     ;
    delta    = delta(notCrit)   ;
    tau      = tau(notCrit)     ;
    deltaMod = (delta - 1).^2   ;
    for k = 55:56
        m        = k - 51	;
        p        = k - 54	;

        Theta       = GetTheta     (deltaMod,tau,A(p),betaInv(m));
        Delta       = GetDelta     (deltaMod,Theta,B(p),a(p));
        Psi         = GetPsi       (deltaMod,tau,C(p),D(p));
        Psi_d       = GetPsi_d     (delta,Psi,C(p));
        Deltabi_d   = GetDeltabi_d (delta,deltaMod,Delta,Theta,A(p),B(p),a(p),b(p),betaInv(m));

        % Calculate zeroth derivative Helmholtz free energy component
        Part    = n(k) * Delta.^(b(p)) .* delta .* Psi  ;
        %         [Helm(notCrit),SumErr(notCrit)] = KahanSum(Helm(notCrit),Part,SumErr(notCrit))      ;
        Helm(notCrit) = Helm(notCrit) + Part;
        
        
        % Calculate first derivative Helmholtz free energy component
        Part        =  Delta.^(b(p)) .* (Psi + delta .* Psi_d)  ;
        Part        = n(k)*(Part + Deltabi_d .* delta .* Psi)   ;
        %         [Helm_d(notCrit),SumErr_d(notCrit)] = KahanSum(Helm_d(notCrit),Part,SumErr_d(notCrit));
        Helm_d(notCrit) = Helm_d(notCrit) + Part; 

    end
    
end
