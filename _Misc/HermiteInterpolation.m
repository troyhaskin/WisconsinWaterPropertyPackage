function y = HermiteInterpolation(x,xp,yp,dydxp)
    
    % Columnify the data vectors
    xp    = xp   (:)    ;
    yp    = yp   (:)    ;
    dydxp = dydxp(:)    ;
    
    % Constants
    N   = length(xp);
    Nsq = N^2       ;
    c   = 2:(Nsq-1) ;
    
    % Allocate Vandermonde arrays
    V   = ones(N,Nsq)   ; % Vandermonde array
    V_x = ones(N,Nsq)   ; % derivative of Vandermonde array
    
    % Zero out first column of derivative matrix
    V_x(:,1) = 0        ;
    
    % Compute the Vandermonde array
    for k = 2:Nsq
        V(:,k) = xp.*V(:,k-1);
    end
    
    % Compute the derivative of Vandermonde array
    for k = 3:Nsq
        V_x(:,k) = c(k-2) * V(:,k-1);
    end
    
    % Form the Hermite polynomial coefficent matrix
    H = [V;V_x];
    
    % Form the solution vector
    b = [yp;dydxp];
    
    % Compute the Hermite polynomial weights.
    w = H\b;
    
    % Use Horner's Method to computer the interpolated values
    y = HornersMethod(x,w(end:-1:1),[0,1]);
    
end


