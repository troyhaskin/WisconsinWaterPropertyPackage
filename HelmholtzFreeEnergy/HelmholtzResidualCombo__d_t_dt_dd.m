function [Helm,Helm_d,Helm_t,Helm_dt,Helm_dd] = HelmholtzResidualCombo__d_t_dt_dd(delta,tau)
    
    persistent c d t n alpha beta gamma epsilon A B C D a b
    if isempty(c)
        [c,d,t,n,alpha,beta,gamma,epsilon,A,B,C,D,a,b] = Coefficients_HelmholtzResidual();
    end
    
    deltaMod  = (delta - 1).^2 + eps(delta) ;
    betaInv   = 1 ./ beta                   ;
    Helm      = 0*delta                     ;
    Helm_d    = Helm                        ;
    Helm_t    = Helm                        ;
    Helm_dt   = Helm                        ;
    Helm_dd   = Helm                        ;
    SumErr    = Helm                        ; % Summation correction (Kahan summation)
    SumErr_d  = Helm                        ; % Summation correction (Kahan summation)
    SumErr_t  = Helm                        ; % Summation correction (Kahan summation)
    SumErr_dt = Helm                        ; % Summation correction (Kahan summation)
    SumErr_dd = Helm                        ; % Summation correction (Kahan summation)
    
    for k = 1:7

        % Shared portion
        Part1 = n(k) * delta.^(d(k)-2) .* tau.^(t(k));


        % Calculate zeroth derivative Helmholtz free energy component
        % Part = n(k) * delta.^(d(k)) .* tau.^(t(k))   ;
        [Helm,SumErr] = KahanSum(Helm, delta.^2 .* Part1 ,SumErr)      ;


        % Calculate first derivative Helmholtz free energy component
        % Part    = n(k) * d(k) * delta.^(d(k)-1) .* tau.^(t(k))  ;
        [Helm_d,SumErr_d] = KahanSum(Helm_d, d(k) * delta .* Part1 ,SumErr_d)      ;

        
        % Calculate first derivative Helmholtz free energy component
        Part2 = t(k) * delta ./ tau .* Part1;
        [Helm_dt,SumErr_dt] = KahanSum(Helm_dt, d(k) * Part2 ,SumErr_dt);


        % Calculate first derivative Helmholtz free energy component
        [Helm_t,SumErr_t] = KahanSum(Helm_t, delta .* Part2 ,SumErr_t);

        
        % Calculate second derivative Helmholtz free energy component
        % Part    = n(k) * d(k) * (d(k)-1) * delta.^(d(k)-2) .* tau.^(t(k))  ;
        [Helm_dd,SumErr_dd] = KahanSum(Helm_dd, d(k) * (d(k)-1) * Part1 ,SumErr_dd);

    end
    
    for k = 8:51

        % Helper variables
        delta2c = delta.^(c(k));
        
        % Calculate zeroth derivative Helmholtz free energy component
        Part1    = n(k) * delta.^(d(k)) .* tau.^(t(k)) .* exp(-delta2c)  ;
        [Helm,SumErr] = KahanSum(Helm,Part1,SumErr)                      ;


        % Calculate first derivative Helmholtz free energy component
        Part_d = Part1 .* (d(k) - c(k) * delta2c) ./ delta   ;
        [Helm_d,SumErr_d] = KahanSum(Helm_d,Part_d,SumErr_d);

        
        %   Calculate
        [Helm_t,SumErr_t] = KahanSum(Helm_t,Part1 * t(k)./tau,SumErr_t);

        
        %   Calculate
        [Helm_dt,SumErr_dt] = KahanSum(Helm_dt,Part_d * t(k)./tau,SumErr_dt);


        % Calculate second derivative Helmholtz free energy component
        Part1 = -(d(k) + c(k)*(c(k)-1)* delta2c) .* Part1  ./ delta.^2   + ...
                 (d(k) - c(k)         * delta2c) .* Part_d ./ delta      ;
        [Helm_dd,SumErr_dd] = KahanSum(Helm_dd,Part1,SumErr_dd)          ;

    end
    
    for k = 52:54
        m       = k - 51                                                        ;
        Arg     = -alpha(m)*(delta-epsilon(m)).^2 - beta(m)*(tau-gamma(m)).^2   ;


        % Calculate zeroth derivative Helmholtz free energy component
        Part1    = n(k) * delta.^(d(k)) .* tau.^(t(k)) .* exp(Arg)   ;
        [Helm,SumErr] = KahanSum(Helm,Part1,SumErr)                  ;


        % Calculate first derivative Helmholtz free energy component
        Part_d = Part1 .* (d(k)./delta - 2*alpha(m)*(delta-epsilon(m)))  ;
        [Helm_d,SumErr_d] = KahanSum(Helm_d,Part_d,SumErr_d)            ;


        %   Calculate
        Part_t = (t(k)./tau - 2.*beta(m).*(tau-gamma(m)));
        [Helm_t,SumErr_t] = KahanSum(Helm_t, Part1 .* Part_t ,SumErr_t);


        %   Calculate
        [Helm_dt,SumErr_dt] = KahanSum(Helm_dt, Part_d .* Part_t ,SumErr_dt);


        % Calculate second derivative Helmholtz free energy component
        Part1 = -(d(k)./delta.^2 + 2*alpha(m))                    .* Part1    + ...
                (d(k)./delta    - 2*alpha(m)*(delta-epsilon(m))) .* Part_d  ;
        [Helm_dd,SumErr_dd] = KahanSum(Helm_dd,Part1,SumErr_dd)              ;
        
    end
    
    for k = 55:56
        m        = k - 51	;
        p        = k - 54	;

        Theta       = GetTheta     (deltaMod,tau,A(p),betaInv(m));
        Delta       = GetDelta     (deltaMod,Theta,B(p),a(p));
        Psi         = GetPsi       (deltaMod,tau,C(p),D(p));
        Psi_d       = GetPsi_d     (delta,Psi,C(p));
        Psi_t       = GetPsi_t     (tau,Psi,D(p));
        Psi_dt      = GetPsi_dt    (delta,tau,Psi,C(p),D(p))    ;
        Psi_dd      = GetPsi_dd    (deltaMod,Psi,C(p));
        Deltabi_d   = GetDeltabi_d (delta,deltaMod,Delta,Theta,A(p),B(p),a(p),b(p),betaInv(m));
        Deltabi_t   = GetDeltabi_t (Delta,Theta,b(p));
        Deltabi_dd  = GetDeltabi_dd(delta,deltaMod,Delta,Theta,A(p),B(p),a(p),b(p),betaInv(m));
        Deltabi_dt  = GetDeltabi_dt(delta,deltaMod,Delta,Theta,A(p),B(p),a(p),b(p),betaInv(m));


        Delta2b = Delta.^(b(p));
        
        % Calculate zeroth derivative Helmholtz free energy component
        Part1    = n(k) * Delta2b .* delta .* Psi  ;
        [Helm,SumErr] = KahanSum(Helm,Part1,SumErr)      ;


        % Calculate first derivative Helmholtz free energy component
        Part1        =  Delta2b .* (Psi + delta .* Psi_d)  ;
        Part1        = n(k)*(Part1 + Deltabi_d .* delta .* Psi)   ;
        [Helm_d,SumErr_d] = KahanSum(Helm_d,Part1,SumErr_d)      ;    
        
        
        %   Calculate
        Part1 = n(k) * delta .*(Deltabi_t .* Psi + Delta2b .* Psi_t);
        [Helm_t,SumErr_t] = KahanSum(Helm_t,Part1,SumErr_t);
        
        
        %   Calculate
        Part1        = Delta2b .* (Psi_t + delta .* Psi_dt)         ;
        Part1        = Part1 + delta .* Deltabi_d .* Psi_t          ;
        Part1        = Part1 + Deltabi_t .*(Psi + delta .* Psi_d)   ;
        Part1        = n(k)*(Part1 + Deltabi_dt .* delta .* Psi)    ;  
        [Helm_dt,SumErr_dt] = KahanSum(Helm_dt,Part1,SumErr_dt)     ;


        % Calculate second derivative Helmholtz free energy component
        Part1      = Delta.^(b(p)).*(2*Psi_d + delta.*Psi_dd)                        ;
        Part1      = Part1 + 2*Deltabi_d .* (Psi+delta.*Psi_d)+Deltabi_dd.*delta.*Psi ;
        [Helm_dd,SumErr_dd] = KahanSum(Helm_dd,Part1,SumErr_dd)             ;

    end
    
end
