function [Helm,Helm_d,Helm_t,Helm_dt,Helm_dd,Helm_tt] = HelmholtzResidualCombo__d_t_dt_dd_tt(delta,tau)
    
    [c,d,t,n,alpha,beta,gamma,epsilon,A,B,C,D,a,b] = Coefficients_HelmholtzResidual();
    
    betaInv = 1 ./ beta     ;
    N       = numel(delta)  ;
    Helm    = zeros(N,1)    ;
    Helm_d  = zeros(N,1)    ;
    Helm_t  = zeros(N,1)    ;
    Helm_dt = zeros(N,1)    ;
    Helm_dd = zeros(N,1)    ;
    Helm_tt = zeros(N,1)    ;
    
    for k = 1:7
        
        Helm0   = n(k) * delta.^d(k) .* tau.^t(k)   ;
        Helm    = Helm + Helm0                      ;
        work    = d(k) * Helm0./delta               ;
        Helm_d  = Helm_d  + work                    ;
        Helm_dd = Helm_dd + (d(k)-1) * work./delta  ;
        work    = t(k) * Helm0./tau                 ;
        Helm_t  = Helm_t  + work                    ;
        Helm_tt = Helm_tt + (t(k)-1) * work./tau    ;
        Helm_dt = Helm_dt + d(k)     * work./delta  ;
        
    end
    
    %     for k = 8:51
    %
    %         % Helper variables
    %         delta2c = delta.^(c(k));
    %
    %         % Calculate zeroth derivative Helmholtz free energy component
    %         Part1 = n(k) * delta.^(d(k)) .* tau.^(t(k)) .* exp(-delta2c);
    %         Helm  = Helm + Part1                                        ;
    %
    %         % Calculate first derivative Helmholtz free energy component
    %         Part_d = Part1 .* (d(k) - c(k) * delta2c) ./ delta  ;
    %         Helm_d = Helm_d + Part_d                            ;
    %
    %         %   Calculate
    %         Helm_t = Helm_t + t(k) * Part1./tau;
    %
    %         %   Calculate
    %         Helm_dt = Helm_dt + t(k) * Part_d./tau;
    %
    %         % Calculate second derivative Helmholtz free energy component
    %         Part1 = -(d(k) + c(k)*(c(k)-1)* delta2c) .* Part1  ./ delta.^2   + ...
    %                  (d(k) - c(k)         * delta2c) .* Part_d ./ delta      ;
    %         Helm_dd = Helm_dd + Part1;
    %
    %     end
    
    
    dd = d(8:51);
    cc = c(8:51);
    tt = t(8:51);
    nn = n(8:51);
    % Helper variables
    delta2c = bsxfun(@power,delta,cc);
    
    % Calculate zeroth derivative Helmholtz free energy component
    Part1 = bsxfun(@times,nn,bsxfun(@power,delta,dd) .* bsxfun(@power,tau,tt) .* exp(-delta2c));
    Helm  = Helm + sum(Part1,2);
    Part_d = Part1.*bsxfun(@rdivide,bsxfun(@minus,dd,bsxfun(@times,cc,delta2c)),delta);
    Helm_d = Helm_d + sum(Part_d,2);
    Helm_t = Helm_t + sum(bsxfun(@times,tt,bsxfun(@rdivide,Part1,tau)),2);
    Helm_dt = Helm_dt + sum(bsxfun(@times,tt,bsxfun(@rdivide,Part_d,tau)),2);
    Helm_dd = Helm_dd + sum(...
        bsxfun(@rdivide,bsxfun(@minus,dd,bsxfun(@times,cc        ,delta2c)).*Part_d,delta) -...
        bsxfun(@rdivide,bsxfun(@plus ,dd,bsxfun(@times,cc.*(cc-1),delta2c)).*Part1 ,delta.^2),2);
    
    
    
    
    
    for k = 52:54
        m       = k - 51                                                        ;
        Arg     = -alpha(m)*(delta-epsilon(m)).^2 - beta(m)*(tau-gamma(m)).^2   ;
        
        
        % Calculate zeroth derivative Helmholtz free energy component
        Part1 = n(k) * delta.^(d(k)) .* tau.^(t(k)) .* exp(Arg) ;
        Helm  = Helm + Part1                                    ;
        %         [Helm,SumErr] = KahanSum(Helm,Part1,SumErr)                  ;
        
        
        % Calculate first derivative Helmholtz free energy component
        Part_d = Part1 .* (d(k)./delta - 2*alpha(m)*(delta-epsilon(m)))  ;
        Helm_d = Helm_d + Part_d;
        %         [Helm_d,SumErr_d] = KahanSum(Helm_d,Part_d,SumErr_d)            ;
        
        
        %   Calculate
        Part_t = (t(k)./tau - 2.*beta(m).*(tau-gamma(m)));
        Helm_t = Helm_t + Part1 .* Part_t;
        %         [Helm_t,SumErr_t] = KahanSum(Helm_t, Part1 .* Part_t ,SumErr_t);
        
        
        %   Calculate
        Helm_dt = Helm_dt + Part_d .* Part_t;
        %         [Helm_dt,SumErr_dt] = KahanSum(Helm_dt, Part_d .* Part_t ,SumErr_dt);
        
        
        % Calculate second derivative Helmholtz free energy component
        Part1 = -(d(k)./delta.^2 + 2*alpha(m))                    .* Part1    + ...
            (d(k)./delta    - 2*alpha(m)*(delta-epsilon(m))) .* Part_d  ;
        Helm_dd = Helm_dd + Part1;
        %         [Helm_dd,SumErr_dd] = KahanSum(Helm_dd,Part1,SumErr_dd)              ;
        
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
        Psi_t       = GetPsi_t     (tau,Psi,D(p));
        Psi_dt      = GetPsi_dt    (delta,tau,Psi,C(p),D(p))    ;
        Psi_dd      = GetPsi_dd    (deltaMod,Psi,C(p));
        [Deltabi_d,Deltabi_t,Deltabi_dt,Deltabi_dd] = ...
            GetDeltabiCombo_d_t_dt_dd(delta,deltaMod,Delta,Theta,A(p),B(p),a(p),b(p),betaInv(m));
        
        
        Delta2b = Delta.^(b(p));
        
        % Calculate zeroth derivative Helmholtz free energy component
        Part1         = n(k) * Delta2b .* delta .* Psi  ;
        Helm(notCrit) = Helm(notCrit) + Part1           ;
        %         [Helm,SumErr] = KahanSum(Helm,Part1,SumErr)      ;
        
        
        % Calculate first derivative Helmholtz free energy component
        Part1           =  Delta2b .* (Psi + delta .* Psi_d)         ;
        Part1           = n(k)*(Part1 + Deltabi_d .* delta .* Psi)   ;
        Helm_d(notCrit) = Helm_d(notCrit) + Part1                    ;
        
        
        %   Calculate
        Part1           = n(k) * delta .*(Deltabi_t .* Psi + Delta2b .* Psi_t)  ;
        Helm_t(notCrit) = Helm_t(notCrit) + Part1                               ;
        %         [Helm_t,SumErr_t] = KahanSum(Helm_t,Part1,SumErr_t);
        
        
        %   Calculate
        Part1            = Delta2b .* (Psi_t + delta .* Psi_dt)      ;
        Part1            = Part1 + delta .* Deltabi_d .* Psi_t       ;
        Part1            = Part1 + Deltabi_t .*(Psi + delta .* Psi_d);
        Part1            = n(k)*(Part1 + Deltabi_dt .* delta .* Psi) ;
        Helm_dt(notCrit) = Helm_dt(notCrit) + Part1                  ;
        %         [Helm_dt,SumErr_dt] = KahanSum(Helm_dt,Part1,SumErr_dt)     ;
        
        
        % Calculate second derivative Helmholtz free energy component
        Part1            = Delta.^(b(p)).*(2*Psi_d + delta.*Psi_dd)                         ;
        Part1            = Part1 + 2*Deltabi_d .* (Psi+delta.*Psi_d)+Deltabi_dd.*delta.*Psi ;
        Helm_dd(notCrit) = Helm_dd(notCrit) + Part1                                         ;
        %         [Helm_dd,SumErr_dd] = KahanSum(Helm_dd,Part1,SumErr_dd)             ;
        
    end
    
end
